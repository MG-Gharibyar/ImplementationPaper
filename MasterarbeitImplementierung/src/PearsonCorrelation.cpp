#include "openfhe.h"
#include <json/json.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// header files needed for serialization
#include "ciphertext-ser.h"
#include "cryptocontext-ser.h"
#include "key/key-ser.h"
#include "scheme/ckksrns/ckksrns-ser.h"

using namespace lbcrypto;

struct FleetData {
    std::vector<double> DriverID;
    std::vector<double> LocLat;
    std::vector<double> LocLong;
    std::vector<double> DistanceTravelled;
    std::vector<double> CO2Emissions;
    std::vector<double> EmissionFactor;
    std::vector<double> Speed;
    std::vector<double> FuelConsumption;
    std::vector<double> CargoWeight;
};

struct FleetDataCipher {
    Ciphertext<DCRTPoly> DriverID;
    Ciphertext<DCRTPoly> LocLat;
    Ciphertext<DCRTPoly> LocLong;
    Ciphertext<DCRTPoly> DistanceTravelled;
    Ciphertext<DCRTPoly> CO2Emissions;
    Ciphertext<DCRTPoly> EmissionFactor;
    Ciphertext<DCRTPoly> Speed;
    Ciphertext<DCRTPoly> FuelConsumption;
    Ciphertext<DCRTPoly> CargoWeight;
};

std::map<int, FleetData> FleetDataMap;

struct Results {
    int batchSize = -1;
    ScalingTechnique RescaleTech = ScalingTechnique::INVALID_RS_TECHNIQUE;
    SecurityLevel SecLevel = SecurityLevel::HEStd_NotSet;
    int FirstModSize = -1;
    int ScaleModSize = -1;
    int RingDimension = -1;
    int MultiplicativeDepth = -1;
    int LevelBudget = -1;
    double Log2q = -1.0;

    double DurationEncryption;
    double DurationCalculation;
    double DurationCalculationPearson;
    double DurationDecryption;
    int LevelForPearson;

    std::string Comment = "";

    double SizeOfCiphertext;
    double ThroughputEncryption;
    double ThroughputDecryption;
    double ThroughputPearsonCalculation;

    int SizeCryptoContext;
    int SizeEvalMultKey;
    int SizeRotKey;

    std::string ActualPearsonResult;
    std::string ExpectedPearsonResult;
    double DurationPearsonCalculation;

    int PrecisionAfterPearsonCalculation;
};

class FHEBenchmark {
private:
    CryptoContext<DCRTPoly> cryptoContext;
    KeyPair<DCRTPoly> keyPair;
    std::map<int, FleetDataCipher> CipherInput;
    std::vector<Ciphertext<DCRTPoly>> PearsonResult;
    uint32_t numSlots;
    usint depth;
    std::vector<double> RawPearson;
    Results BenchmarkResults;

    std::vector<double> CalcPlaintextPearson() {
        int n = FleetDataMap.begin()->second.CargoWeight.size();

        std::vector<double> result;

        for (const auto& [truckId, data] : FleetDataMap) {
            double meanFuelConsumption = 0.0, meanCargoWeight = 0.0;
            for (int i = 0; i < n; ++i) {
                meanFuelConsumption += data.FuelConsumption[i];
                meanCargoWeight += data.CargoWeight[i];
            }

            meanFuelConsumption /= n;
            meanCargoWeight /= n;

            double sum_product_deviations = 0.0, sum_squared_fuelcons = 0.0, sum_squared_cargow = 0.0;
            for (int i = 0; i < n; ++i) {
                double devFuelCons = data.FuelConsumption[i] - meanFuelConsumption;
                double devCargoW = data.CargoWeight[i] - meanCargoWeight;
                sum_product_deviations += devFuelCons * devCargoW;
                sum_squared_fuelcons += devFuelCons * devFuelCons;
                sum_squared_cargow += devCargoW * devCargoW;
            }

            double correlation_coefficient = sum_product_deviations / (sqrt(sum_squared_fuelcons) * sqrt(sum_squared_cargow));
            result.push_back(correlation_coefficient);
        }

        return result;
    }

    double CalculateApproximationError(const std::vector<double>& result, const std::vector<double>& expectedResult) {
        if (result.size() != expectedResult.size())
            OPENFHE_THROW(config_error, "Cannot compare vectors with different numbers of elements");

        double maxError = 0;
        for (size_t i = 0; i < result.size(); ++i) {
            double error = std::abs(result[i] - expectedResult[i]);
            if (maxError < error)
                maxError = error;
        }

        return std::log2(maxError);
    }

    void printVector(std::string FirstText, const std::vector<double>& vec, int numElements = -1) {
        std::cout << FirstText;
        int count = 0;
        for (const auto& element : vec) {
            if (numElements == -1 || count < numElements) {
                std::cout << std::fixed << std::setprecision(7) << element << " ";
                count++;
            } else {
                break;
            }
        }
        std::cout << std::endl;
    }

    std::string vectorToFormattedString(const std::vector<double>& vec, int length) {
        std::ostringstream oss;
        for (int i = 0; i < length; ++i) {
            oss << std::fixed << std::setprecision(10) << vec[i] << " ";
        }
        return oss.str();
    }

public:
    FHEBenchmark(ScalingTechnique RescaleTech, SecurityLevel SecLevel, int batchSize, int FirstModsize, int ScaleModSize) {
        CCParams<CryptoContextCKKSRNS> parameters;
        SecretKeyDist secretKeyDist = UNIFORM_TERNARY;

        parameters.SetSecretKeyDist(secretKeyDist);
        parameters.SetSecurityLevel(SecLevel);
        if (SecLevel == SecurityLevel::HEStd_NotSet) { 
            parameters.SetRingDim(1024);
        }

        depth = 12;
        parameters.SetBatchSize(batchSize);
        parameters.SetMultiplicativeDepth(depth);
        parameters.SetScalingModSize(ScaleModSize);
        parameters.SetFirstModSize(FirstModsize);

        cryptoContext = GenCryptoContext(parameters);

        cryptoContext->Enable(PKE);
        cryptoContext->Enable(KEYSWITCH);
        cryptoContext->Enable(LEVELEDSHE);
        cryptoContext->Enable(ADVANCEDSHE);
        cryptoContext->Enable(FHE);

        numSlots = batchSize;

        keyPair = cryptoContext->KeyGen();
        cryptoContext->EvalMultKeyGen(keyPair.secretKey);
        cryptoContext->EvalSumKeyGen(keyPair.secretKey, keyPair.publicKey);
        cryptoContext->EvalRotateKeyGen(keyPair.secretKey, {1, 2, 3, -1, -2, -3});

        RawPearson = CalcPlaintextPearson();

        BenchmarkResults.batchSize = batchSize;
        BenchmarkResults.RescaleTech = RescaleTech;
        BenchmarkResults.SecLevel = SecLevel;
        BenchmarkResults.FirstModSize = FirstModsize;
        BenchmarkResults.ScaleModSize = ScaleModSize;
        BenchmarkResults.RingDimension = cryptoContext->GetRingDimension();
        BenchmarkResults.MultiplicativeDepth = depth;
        double log2q = log2(cryptoContext->GetCryptoParameters()->GetElementParams()->GetModulus().ConvertToDouble());
        BenchmarkResults.Log2q = std::isinf(log2q) ? 10000.0 : log2q;

        std::stringstream s1, s2, s3;
        size_t lengthInBytes;
        Serial::Serialize(cryptoContext, s1, SerType::BINARY);
        std::string str = s1.str();
        BenchmarkResults.SizeCryptoContext = sizeof(str[0]) * (str.length() + 1);

        cryptoContext->SerializeEvalMultKey(s2, SerType::BINARY);
        str = s2.str();
        lengthInBytes = sizeof(str[0]) * (str.length() + 1);
        BenchmarkResults.SizeEvalMultKey = lengthInBytes / 1024;

        cryptoContext->SerializeEvalAutomorphismKey(s3, SerType::BINARY);
        str = s3.str();
        lengthInBytes = sizeof(str[0]) * (str.length() + 1);
        BenchmarkResults.SizeRotKey = lengthInBytes / 1024;
    }

    void EncryptData() {
        TimeVar tStart, t;
        std::cout << "Start of encryption" << std::endl;
        double durEncryption = 0.0;
        TIC(tStart);

        for (const auto& [truckId, data] : FleetDataMap) {
            TIC(t);
            FleetDataCipher& datacipher = CipherInput[truckId];

            Plaintext ptxt = cryptoContext->MakeCKKSPackedPlaintext(data.CargoWeight);
            ptxt->SetLength(numSlots);
            TIC(t);
            datacipher.CargoWeight = cryptoContext->Encrypt(keyPair.publicKey, ptxt);
            durEncryption += TOC_MS(t);

            TIC(t);
            ptxt = cryptoContext->MakeCKKSPackedPlaintext(data.CO2Emissions);
            ptxt->SetLength(numSlots);
            datacipher.CO2Emissions = cryptoContext->Encrypt(keyPair.publicKey, ptxt);
            durEncryption += TOC_MS(t);

            TIC(t);
            ptxt = cryptoContext->MakeCKKSPackedPlaintext(data.DistanceTravelled);
            ptxt->SetLength(numSlots);
            datacipher.DistanceTravelled = cryptoContext->Encrypt(keyPair.publicKey, ptxt);
            durEncryption += TOC_MS(t);

            TIC(t);
            ptxt = cryptoContext->MakeCKKSPackedPlaintext(data.DriverID);
            ptxt->SetLength(numSlots);
            datacipher.DriverID = cryptoContext->Encrypt(keyPair.publicKey, ptxt);
            durEncryption += TOC_MS(t);

            TIC(t);
            ptxt = cryptoContext->MakeCKKSPackedPlaintext(data.EmissionFactor);
            ptxt->SetLength(numSlots);
            datacipher.EmissionFactor = cryptoContext->Encrypt(keyPair.publicKey, ptxt);
            durEncryption += TOC_MS(t);

            TIC(t);
            ptxt = cryptoContext->MakeCKKSPackedPlaintext(data.FuelConsumption);
            ptxt->SetLength(numSlots);
            datacipher.FuelConsumption = cryptoContext->Encrypt(keyPair.publicKey, ptxt);
            durEncryption += TOC_MS(t);

            TIC(t);
            ptxt = cryptoContext->MakeCKKSPackedPlaintext(data.LocLat);
            ptxt->SetLength(numSlots);
            datacipher.LocLat = cryptoContext->Encrypt(keyPair.publicKey, ptxt);
            durEncryption += TOC_MS(t);

            TIC(t);
            ptxt = cryptoContext->MakeCKKSPackedPlaintext(data.LocLong);
            ptxt->SetLength(numSlots);
            datacipher.LocLong = cryptoContext->Encrypt(keyPair.publicKey, ptxt);
            durEncryption += TOC_MS(t);

            TIC(t);
            ptxt = cryptoContext->MakeCKKSPackedPlaintext(data.Speed);
            ptxt->SetLength(numSlots);
            datacipher.Speed = cryptoContext->Encrypt(keyPair.publicKey, ptxt);
            durEncryption += TOC_MS(t);
        }
        double totalDurationMS = TOC_MS(tStart);

        std::cout << "Size of FleetDataMap: " << FleetDataMap.size() << std::endl;
        BenchmarkResults.DurationEncryption = totalDurationMS;

        std::cout << "Total encryption duration: " << totalDurationMS << " ms" << std::endl;
        std::cout << "Average encryption duration per data field: " << (durEncryption / (9 * FleetDataMap.size())) << " ms" << std::endl;
    }

    void SetDataSize() {
        std::stringstream s;
        Serial::Serialize(CipherInput[0].FuelConsumption, s, SerType::BINARY);
        std::string str = s.str();
        size_t lengthInBytes = sizeof(str[0]) * (str.length() + 1);
        BenchmarkResults.SizeOfCiphertext = lengthInBytes / 1024;
    }

    void CalcCipherPearson() {
        TimeVar t, tmean, trsq, tdiv, tsq, tdecrypt;
        TIC(t);
        std::cout << "\nStart of Pearson correlation coefficients calculation" << std::endl;
        std::cout << "Level: " << CipherInput.begin()->second.FuelConsumption->GetLevel() << std::endl;

        std::vector<Ciphertext<DCRTPoly>> pcorresults;

        BenchmarkResults.LevelForPearson = CipherInput.begin()->second.FuelConsumption->GetLevel();

        for (const auto& [truckId, data] : CipherInput) {
            TIC(tmean);
            auto csumfuel = cryptoContext->EvalSum(data.FuelConsumption, numSlots);
            auto csumweight = cryptoContext->EvalSum(data.CargoWeight, numSlots);
            auto cmeanfuel = cryptoContext->EvalMult(0.03125, csumfuel);
            auto cmeanweight = cryptoContext->EvalMult(0.03125, csumweight);
            BenchmarkResults.DurationCalculation += RoundToDecimals(TOC_MS(tmean), 2);

            auto cxixavg = cryptoContext->EvalSub(data.FuelConsumption, cmeanfuel);
            auto cyiyavg = cryptoContext->EvalSub(data.CargoWeight, cmeanweight);
            auto numerator = cryptoContext->EvalSum(cryptoContext->EvalMult(cxixavg, cyiyavg), numSlots);

            TIC(tsq);
            auto denom1 = cryptoContext->EvalSum(cryptoContext->EvalSquare(cxixavg), numSlots);
            auto denom2 = cryptoContext->EvalSum(cryptoContext->EvalSquare(cyiyavg), numSlots);
            auto denom3 = cryptoContext->EvalMult(denom1, denom2);
            BenchmarkResults.DurationCalculation += RoundToDecimals(TOC_MS(tsq), 2);

            TIC(trsq);
            auto denom4 = cryptoContext->EvalChebyshevFunction([](double x) -> double { return std::sqrt(x); }, denom3, 5e6, 16e6, 5);
            BenchmarkResults.DurationCalculation += RoundToDecimals(TOC_MS(trsq), 2);

            TIC(tdiv);
            auto denominator = cryptoContext->EvalDivide(denom4, std::sqrt(5e6), std::sqrt(16e6), 5);
            auto result = cryptoContext->EvalMult(numerator, denominator);
            BenchmarkResults.DurationCalculation += RoundToDecimals(TOC_MS(tdiv), 2);

            pcorresults.push_back(result);
        }

        BenchmarkResults.DurationCalculationPearson = RoundToDecimals(TOC_MS(t), 2);
        BenchmarkResults.LevelForPearson = pcorresults[0]->GetLevel();

        std::cout << "Level after Pearson Correlation: " << pcorresults[0]->GetLevel() << std::endl;
        std::cout << "Total calculation duration: " << RoundToDecimals(BenchmarkResults.DurationCalculationPearson, 2) << " ms" << std::endl;

        TIC(tdecrypt);

        try {
            std::vector<double> results;

            for (size_t i = 0; i < pcorresults.size(); i++) {
                Plaintext resultPCor;
                cryptoContext->Decrypt(keyPair.secretKey, pcorresults[i], &resultPCor);
                resultPCor->SetLength(numSlots);
                results.push_back(resultPCor->GetRealPackedValue()[0]);
            }

            BenchmarkResults.DurationDecryption = RoundToDecimals(TOC_MS(tdecrypt), 2);

            BenchmarkResults.PrecisionAfterPearsonCalculation = std::floor(CalculateApproximationError(RawPearson, results));
            BenchmarkResults.ExpectedPearsonResult = vectorToFormattedString(RawPearson, RawPearson.size());
            BenchmarkResults.ActualPearsonResult = vectorToFormattedString(results, results.size());

            std::cout << "\nResults of Pearson correlation" << std::endl;
            printVector("\tExpected result:\t", RawPearson, RawPearson.size());
            printVector("\tActual result:\t\t", results, results.size());
            std::cout << "\tPrecision: " << BenchmarkResults.PrecisionAfterPearsonCalculation << " bits" << std::endl;

        } catch (const std::exception& e) {
            BenchmarkResults.Comment = "Error while decrypting Pearson correlation coefficients results";
            std::cerr << e.what() << '\n';
        }
    }

    double RoundToDecimals(double value, int decimals) {
        if (!std::isfinite(value)) {
            return 0.0;
        }
        double factor = std::pow(10.0, decimals);
        return std::round(value * factor) / factor;
    }

    void CalculateThroughput() {
        int length = CipherInput.size();

        double throughputEncryption = (BenchmarkResults.SizeOfCiphertext * length) / (static_cast<double>(BenchmarkResults.DurationEncryption) / 1000);
        double throughputDecryption = (BenchmarkResults.SizeOfCiphertext * length) / (static_cast<double>(BenchmarkResults.DurationDecryption) / 1000);
        double throughputPearsonCalculation = (BenchmarkResults.SizeOfCiphertext * length) / (static_cast<double>(BenchmarkResults.DurationCalculationPearson) / 1000);

        BenchmarkResults.ThroughputEncryption = RoundToDecimals(throughputEncryption, 2);
        BenchmarkResults.ThroughputDecryption = RoundToDecimals(throughputDecryption, 2);
        BenchmarkResults.ThroughputPearsonCalculation = RoundToDecimals(throughputPearsonCalculation, 2);
    }

    void WriteJSON(const std::string& filename) {
        Json::Value root;
        std::ifstream inputFile(filename);
        if (inputFile.is_open()) {
            inputFile >> root;
            inputFile.close();
        } else {
            root = Json::Value(Json::arrayValue);
            std::cout << "No JSON file found. Creating a new one" << std::endl;
        }

        Json::Value newResult(Json::objectValue);

        newResult["SecurityLevel [bit]"] = static_cast<int>(BenchmarkResults.SecLevel);
        newResult["ScalingTechnique"] = static_cast<int>(BenchmarkResults.RescaleTech);
        newResult["FirstModSize [bit]"] = BenchmarkResults.FirstModSize;
        newResult["ScaleModSize [bit]"] = BenchmarkResults.ScaleModSize;
        newResult["MultiplicativeDepth [L]"] = BenchmarkResults.MultiplicativeDepth;
        newResult["RingDimension [bit]"] = BenchmarkResults.RingDimension;
        newResult["Log2(q) [-]"] = BenchmarkResults.Log2q;
        newResult["BatchSize"] = BenchmarkResults.batchSize;
        newResult["LevelBudget"] = BenchmarkResults.LevelBudget;

        newResult["DurationCalculationPearson [ms]"] = BenchmarkResults.DurationCalculationPearson;
        newResult["PrecisionAfterPearsonCalculation [bit]"] = BenchmarkResults.PrecisionAfterPearsonCalculation;
        newResult["ActualPearsonResult"] = BenchmarkResults.ActualPearsonResult;
        newResult["ExpectedPearsonResult"] = BenchmarkResults.ExpectedPearsonResult;
        newResult["Duration Decryption [ms]"] = BenchmarkResults.DurationDecryption;
        newResult["Duration Encryption [ms]"] = BenchmarkResults.DurationEncryption;

        newResult["SizeCryptoContext [byte]"] = BenchmarkResults.SizeCryptoContext;
        newResult["SizeOfCiphertext [kb]"] = BenchmarkResults.SizeOfCiphertext;
        newResult["SizeEvalMultKey [kb]"] = BenchmarkResults.SizeEvalMultKey;
        newResult["SizeRotKey [kb]"] = BenchmarkResults.SizeRotKey;

        newResult["ThroughputDecryption [kb/s]"] = RoundToDecimals(BenchmarkResults.ThroughputDecryption, 2);
        newResult["ThroughputEncryption [kb/s]"] = RoundToDecimals(BenchmarkResults.ThroughputEncryption, 2);
        newResult["ThroughputPearsonCalculation [kb/s]"] = RoundToDecimals(BenchmarkResults.ThroughputPearsonCalculation, 2);

        newResult["Decryption Comment"] = BenchmarkResults.Comment;

        root.append(newResult);

        std::ofstream outputFile(filename);
        if (outputFile.is_open()) {
            Json::StreamWriterBuilder writer;
            writer["indentation"] = "    ";
            std::unique_ptr<Json::StreamWriter> jsonWriter(writer.newStreamWriter());
            jsonWriter->write(root, &outputFile);
            outputFile.close();
            std::cout << "Results appended to JSON file." << std::endl;
        } else {
            std::cerr << "Error while writing JSON file." << std::endl;
        }
    }

    void ClearFHE() {
        cryptoContext->ClearEvalMultKeys();
        cryptoContext->ClearEvalAutomorphismKeys();
        lbcrypto::CryptoContextFactory<lbcrypto::DCRTPoly>::ReleaseAllContexts();
    }
};

void Benchmark(
    const std::string& filename,
    int repetitions,
    ScalingTechnique startScalingTech,
    SecurityLevel startSecLevel,
    int startBatchSize,
    int startFirstModSize,
    int startScaleModSize,
    bool startFromGivenParams = false
) {
    bool foundStart = !startFromGivenParams; // Start immediately if not using specific parameters
    int i = 0;

    for (ScalingTechnique ST : {ScalingTechnique::FIXEDAUTO, ScalingTechnique::FLEXIBLEAUTO, ScalingTechnique::FLEXIBLEAUTOEXT}) {
        for (SecurityLevel SL : {SecurityLevel::HEStd_128_classic, SecurityLevel::HEStd_192_classic, SecurityLevel::HEStd_256_classic, SecurityLevel::HEStd_128_quantum, SecurityLevel::HEStd_192_quantum, SecurityLevel::HEStd_256_quantum}) {
            for (int BS : {32}) {
                for (auto ModSize : {std::make_pair(60, 30) ,std::make_pair(60, 40), std::make_pair(60, 50), std::make_pair(60, 59)}) {
                    
                    // Check if the current parameters match the start parameters
                    if (!foundStart) {
                        if (ST == startScalingTech && SL == startSecLevel && BS == startBatchSize &&
                            ModSize.first == startFirstModSize && ModSize.second == startScaleModSize) {
                            foundStart = true;  // Found the start point
                        } else {
                            continue;  // Skip until we find the start point
                        }
                    }

                    std::cout << "##############################" << std::endl;
                    std::cout << "  Benchmark no " << ++i << std::endl;
                    std::cout << "##############################" << std::endl;

                    for (int rep = 0; rep < repetitions; rep++) {
                        std::cout << "\n====== Repetition " << rep + 1 << "/" << repetitions << " =====\n";
                        FHEBenchmark bm(ST, SL, BS, ModSize.first, ModSize.second);
                        bm.EncryptData();
                        bm.SetDataSize();
                        bm.CalcCipherPearson();
                        bm.CalculateThroughput();
                        bm.WriteJSON(filename);
                        bm.ClearFHE();
                    }
                }
            }
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc != 8) {
        std::cerr << "Usage: " << argv[0] << " <Filename.json> <Repetitions> <ScalingTechnique> <SecurityLevel> <BatchSize> <FirstModSize> <ScaleModSize>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    int repetitions = std::atoi(argv[2]);
    ScalingTechnique startScalingTech = static_cast<ScalingTechnique>(std::atoi(argv[3]));
    SecurityLevel startSecLevel = static_cast<SecurityLevel>(std::atoi(argv[4]));
    int startBatchSize = std::atoi(argv[5]);
    int startFirstModSize = std::atoi(argv[6]);
    int startScaleModSize = std::atoi(argv[7]);

    if (repetitions <= 0) {
        std::cerr << "Repetitions must be a positive integer" << std::endl;
        return 1;
    }

    Json::Value root;
    Json::CharReaderBuilder builder;
    std::ifstream file("FleetData.json");
    std::string errors;

    if (!Json::parseFromStream(builder, file, &root, &errors)) {
        std::cout << "Error while parsing json file: " << errors << std::endl;
    }

    for (const auto& item : root) {
        int truckID = item["TruckID"].asInt();
        FleetData& data = FleetDataMap[truckID];

        data.DriverID.push_back(item["DriverID"].asDouble());
        data.LocLat.push_back(item["LocLat"].asDouble());
        data.LocLong.push_back(item["LocLong"].asDouble());
        data.DistanceTravelled.push_back(item["DistanceTravelled"].asDouble());
        data.CO2Emissions.push_back(item["CO2Emissions"].asDouble());
        data.EmissionFactor.push_back(item["EmissionFactor"].asDouble());
        data.Speed.push_back(item["Speed"].asDouble());
        data.FuelConsumption.push_back(item["FuelConsumption"].asDouble());
        data.CargoWeight.push_back(item["CargoWeight"].asDouble());
    }

    TimeVar t;
    TIC(t);


    Benchmark(filename, repetitions, startScalingTech, startSecLevel, startBatchSize, startFirstModSize, startScaleModSize, true);

    int milliseconds = TOC(t);
    int seconds = milliseconds / 1000;
    int hours = seconds / 3600;
    int remainingSeconds = seconds % 3600;
    int minutes = remainingSeconds / 60;

    std::cout << "\nBenchmarking completed in " << hours << " hours and " << minutes << " minutes." << std::endl;

    return 0;
}