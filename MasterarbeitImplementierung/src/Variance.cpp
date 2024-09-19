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

struct FleetData
{
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

struct FleetDataCipher
{
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

std::map<int, FleetData> FleetDataMap; // Die Map enthält alle Daten aus der Funktion ReadDataFromJson()

struct Results
{
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
    double DurationCalculationVariance;
    double DurationDecryption;
    int LevelForVariance;

    std::string Comment = "";

    double SizeOfCiphertext;
    double ThroughputEncryption;
    double ThroughputDecryption;
    double ThroughputVarianceCalculation;

    int SizeCryptoContext;
    int SizeEvalMultKey;
    int SizeRotKey;

    std::string ActualVarianceResult;
    std::string ExpectedVarianceResult;
    double DurationVarianceCalculation;

    int PrecisionAfterVarianceCalculation;
};

class FHEBenchmark
{
private:
    CryptoContext<DCRTPoly> cryptoContext;
    KeyPair<DCRTPoly> keyPair;
    std::map<int, FleetDataCipher> CipherInput;
    std::vector<Ciphertext<DCRTPoly>> VarianceResult;
    uint32_t numSlots;
    usint depth;
    std::vector<double> RawVariance;
    Results BenchmarkResults;

    std::vector<double> CalcPlaintextVariance()
    {
        int n = FleetDataMap.begin()->second.CargoWeight.size(); // number of data points
        std::vector<double> result;

        for (const auto &[truckId, data] : FleetDataMap)
        {
            double mean = 0.0;
            for (int i = 0; i < n; i++)
            {
                mean += data.FuelConsumption[i];
            }
            mean /= n;

            // Sum (x_i - µ)^2 calculation
            double sum = 0.0;
            for (int i = 0; i < n; i++)
            {
                sum += std::pow((data.FuelConsumption[i] - mean), 2);
            }

            sum /= n;
            result.push_back(sum);
        }

        return result;
    }

    double CalculateApproximationError(const std::vector<double> &result, const std::vector<double> &expectedResult)
    {
        if (result.size() != expectedResult.size())
            OPENFHE_THROW(config_error, "Cannot compare vectors with different numbers of elements");

        double maxError = 0;
        for (size_t i = 0; i < result.size(); ++i)
        {
            double error = std::abs(result[i] - expectedResult[i]);
            if (maxError < error)
                maxError = error;
        }

        return std::log2(maxError);
    }

    double CalculateApproximationError(const std::vector<double> &result, const std::vector<std::complex<double>> &expectedResult)
    {
        double maxError = 0;
        for (size_t i = 0; i < result.size(); ++i)
        {
            double error = std::abs(result[i] - expectedResult[i].real());
            if (maxError < error)
                maxError = error;
        }

        double prec = std::log2(maxError);
        return std::isfinite(prec) ? prec : -1.0;
    }

    double CalculateApproximationError(const std::vector<double> &result, const std::vector<std::complex<double>> &expectedResult, int length)
    {
        double maxError = 0;
        for (int i = 0; i < length; ++i)
        {
            double error = std::abs(result[i] - expectedResult[i].real());
            if (maxError < error)
                maxError = error;
        }

        double prec = std::log2(maxError);
        return std::isfinite(prec) ? prec : -1.0;
    }

    void printVector(std::string FirstText, const std::vector<double> &vec, int numElements = -1)
    {
        std::cout << FirstText;
        int count = 0;
        for (const auto &element : vec)
        {
            if (numElements == -1 || count < numElements)
            {
                std::cout << std::fixed << std::setprecision(7) << element << " ";
                count++;
            }
            else
            {
                break;
            }
        }
        std::cout << std::endl;
    }

    std::vector<double> stringToDoubleVector(const std::string &str)
    {
        std::vector<double> result;
        std::stringstream ss(str);
        double value;
        while (ss >> value)
        {
            result.push_back(value);
            while (ss.peek() == ' ' || ss.peek() == ',' || ss.peek() == '\t')
            {
                ss.ignore();
            }
        }
        return result;
    }

    std::string vectorToFormattedString(const std::vector<double> &vec, int length)
    {
        std::ostringstream oss;
        for (int i = 0; i < length; ++i)
        {
            oss << std::fixed << std::setprecision(10) << vec[i] << " ";
        }
        return oss.str();
    }

public:
    FHEBenchmark(ScalingTechnique RescaleTech, SecurityLevel SecLevel, int batchSize, int FirstModsize, int ScaleModSize)
    {
        CCParams<CryptoContextCKKSRNS> parameters;
        SecretKeyDist secretKeyDist = UNIFORM_TERNARY;

        parameters.SetSecretKeyDist(secretKeyDist);
        parameters.SetSecurityLevel(SecLevel);
        if (SecLevel == SecurityLevel::HEStd_NotSet)
        {
            parameters.SetRingDim(1024);
        }

        parameters.SetBatchSize(batchSize);
        parameters.SetMultiplicativeDepth(3);
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

        RawVariance = CalcPlaintextVariance();

        BenchmarkResults.batchSize = batchSize;
        BenchmarkResults.RescaleTech = RescaleTech;
        BenchmarkResults.SecLevel = SecLevel;
        BenchmarkResults.FirstModSize = FirstModsize;
        BenchmarkResults.ScaleModSize = ScaleModSize;
        BenchmarkResults.RingDimension = cryptoContext->GetRingDimension();
        BenchmarkResults.MultiplicativeDepth = 3;
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

    void EncryptData()
    {
        TimeVar tStart, t;
        std::cout << "Start of encryption" << std::endl;
        double durEncryption = 0.0;
        TIC(tStart);

        for (const auto &[truckId, data] : FleetDataMap)
        {
            TIC(t);
            FleetDataCipher &datacipher = CipherInput[truckId];

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

    void SetDataSize()
    {
        std::stringstream s;
        Serial::Serialize(CipherInput[0].FuelConsumption, s, SerType::BINARY);
        std::string str = s.str();
        size_t lengthInBytes = sizeof(str[0]) * (str.length() + 1);
        BenchmarkResults.SizeOfCiphertext = lengthInBytes / 1024;
    }

    void CalcCipherVariance()
    {
        TimeVar t;
        TIC(t);
        std::cout << "\nStart of variance calculation" << std::endl;

        std::vector<Ciphertext<DCRTPoly>> varresults;
        BenchmarkResults.LevelForVariance = CipherInput.begin()->second.FuelConsumption->GetLevel();
        std::cout << "\tLevel before variance calculation: " << CipherInput.begin()->second.FuelConsumption->GetLevel() << std::endl;

        for (const auto &[truckId, data] : CipherInput)
        {
            auto csum = cryptoContext->EvalSum(data.FuelConsumption, numSlots);
            auto cmean = cryptoContext->EvalMult(0.03125, csum);

            auto cdiff = cryptoContext->EvalSub(data.FuelConsumption, cmean);
            auto csq = cryptoContext->EvalSquare(cdiff);

            auto csum2 = cryptoContext->EvalSum(csq, numSlots);
            auto cresult = cryptoContext->EvalMult(0.03125, csum2);

            varresults.push_back(cresult);
        }

        BenchmarkResults.DurationCalculationVariance = TOC_MS(t);
        BenchmarkResults.LevelForVariance = varresults[0]->GetLevel();

        std::cout << "Level after variance calculation: " << varresults[0]->GetLevel() << std::endl;

        try
        {
            std::vector<double> results;
            double totalDecryptionTime = 0.0;

            for (size_t i = 0; i < varresults.size(); i++)
            {
                TimeVar tdecrypt;
                TIC(tdecrypt);

                Plaintext resultVariance;
                cryptoContext->Decrypt(keyPair.secretKey, varresults[i], &resultVariance);

                double decryptionTime = TOC_MS(tdecrypt);
                totalDecryptionTime += decryptionTime;

                resultVariance->SetLength(numSlots);
                results.push_back(resultVariance->GetRealPackedValue()[0]);
            }

            BenchmarkResults.DurationDecryption = totalDecryptionTime;

            double averageDecryptionTime = totalDecryptionTime / varresults.size();
            BenchmarkResults.PrecisionAfterVarianceCalculation = std::floor(CalculateApproximationError(RawVariance, results));
            BenchmarkResults.ExpectedVarianceResult = vectorToFormattedString(RawVariance, RawVariance.size());
            BenchmarkResults.ActualVarianceResult = vectorToFormattedString(results, results.size());

            std::cout << "\nResults of variance calculation" << std::endl;
            printVector("\tExpected result:\t", RawVariance, RawVariance.size());
            printVector("\tActual result:\t\t", results, results.size());
            std::cout << "\tPrecision: " << BenchmarkResults.PrecisionAfterVarianceCalculation << " bits" << std::endl;
            std::cout << "\tTotal decryption duration: " << BenchmarkResults.DurationDecryption << " ms" << std::endl;
            std::cout << "\tAverage decryption time per element: " << averageDecryptionTime << " ms" << std::endl;
        }
        catch (const std::exception &e)
        {
            BenchmarkResults.Comment = "Error while decrypting variance results";
            std::cerr << e.what() << '\n';
        }

        std::cout << "Total variance calculation duration: " << BenchmarkResults.DurationCalculationVariance << " ms" << std::endl;
    }

    double RoundToDecimals(double value, int decimals)
    {
        if (!std::isfinite(value))
        {
            return 0.0;
        }
        double factor = std::pow(10.0, decimals);
        return std::round(value * factor) / factor;
    }

    void CalculateThroughput()
    {
        int length = CipherInput.size();

        double throughputEncryption = (BenchmarkResults.SizeOfCiphertext * length) / (static_cast<double>(BenchmarkResults.DurationEncryption) / 1000);
        double throughputDecryption = (BenchmarkResults.SizeOfCiphertext * length) / (static_cast<double>(BenchmarkResults.DurationDecryption) / 1000);
        double throughputVarianceCalculation = (BenchmarkResults.SizeOfCiphertext * length) / (static_cast<double>(BenchmarkResults.DurationCalculationVariance) / 1000);

        BenchmarkResults.ThroughputEncryption = RoundToDecimals(throughputEncryption, 2);
        BenchmarkResults.ThroughputDecryption = RoundToDecimals(throughputDecryption, 2);
        BenchmarkResults.ThroughputVarianceCalculation = RoundToDecimals(throughputVarianceCalculation, 2);
    }

    void WriteJSON(const std::string &filename)
    {
        Json::Value root;
        std::ifstream inputFile(filename);
        if (inputFile.is_open())
        {
            inputFile >> root;
            inputFile.close();
        }
        else
        {
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

        newResult["DurationCalculationVariance [ms]"] = BenchmarkResults.DurationCalculationVariance;
        newResult["PrecisionAfterVarianceCalculation [bit]"] = BenchmarkResults.PrecisionAfterVarianceCalculation;
        newResult["ActualVarianceResult"] = BenchmarkResults.ActualVarianceResult;
        newResult["ExpectedVarianceResult"] = BenchmarkResults.ExpectedVarianceResult;
        newResult["Duration Decryption [ms]"] = BenchmarkResults.DurationDecryption;
        newResult["Duration Encryption [ms]"] = BenchmarkResults.DurationEncryption;

        newResult["SizeCryptoContext [byte]"] = BenchmarkResults.SizeCryptoContext;
        newResult["SizeOfCiphertext [kb]"] = BenchmarkResults.SizeOfCiphertext;
        newResult["SizeEvalMultKey [kb]"] = BenchmarkResults.SizeEvalMultKey;
        newResult["SizeRotKey [kb]"] = BenchmarkResults.SizeRotKey;

        newResult["ThroughputDecryption [kb/s]"] = RoundToDecimals(BenchmarkResults.ThroughputDecryption, 2);
        newResult["ThroughputEncryption [kb/s]"] = RoundToDecimals(BenchmarkResults.ThroughputEncryption, 2);
        newResult["ThroughputVarianceCalculation [kb/s]"] = RoundToDecimals(BenchmarkResults.ThroughputVarianceCalculation, 2);

        newResult["Decryption Comment"] = BenchmarkResults.Comment;

        root.append(newResult);

        std::ofstream outputFile(filename);
        if (outputFile.is_open())
        {
            Json::StreamWriterBuilder writer;
            writer["indentation"] = "    "; // Optional: Einrückung für lesbare Ausgabe
            std::unique_ptr<Json::StreamWriter> jsonWriter(writer.newStreamWriter());
            jsonWriter->write(root, &outputFile);
            outputFile.close();
            std::cout << "Results appended to JSON file." << std::endl;
        }
        else
        {
            std::cerr << "Error while writing JSON file." << std::endl;
        }
    }

    void ClearFHE()
    {
        cryptoContext->ClearEvalMultKeys();
        cryptoContext->ClearEvalAutomorphismKeys();
        lbcrypto::CryptoContextFactory<lbcrypto::DCRTPoly>::ReleaseAllContexts();
    }
};

void Benchmark(const std::string &filename, int repetitions, ScalingTechnique RescaleTech, SecurityLevel SecLevel, int batchSize, int FirstModsize, int ScaleModSize)
{
    for (int i = 0; i < repetitions; i++)
    {
        std::cout << "\n====== Benchmark (" << RescaleTech << "/" << SecLevel << "/" << batchSize << "/" << FirstModsize << "/" << ScaleModSize << ") =====\n";
        std::cout << "       Repetition " << i + 1 << "/" << repetitions << std::endl;

        FHEBenchmark bm(RescaleTech, SecLevel, batchSize, FirstModsize, ScaleModSize);
        bm.EncryptData();
        bm.SetDataSize();
        bm.CalcCipherVariance();
        // bm.Decrypt();
        bm.CalculateThroughput();
        bm.WriteJSON(filename);
        bm.ClearFHE();
    }
}

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " <Filename.json> <Repetitions>" << std::endl;
        return 1;
    }

    std::string Filename = argv[1];
    int Repetitions = std::atoi(argv[2]);

    if (Repetitions <= 0)
    {
        std::cerr << "Reptitions must be a positive integer" << std::endl;
        return 1;
    }

    Json::Value root;
    Json::CharReaderBuilder builder;
    std::ifstream file("FleetData.json");
    std::string errors;

    if (!Json::parseFromStream(builder, file, &root, &errors))
    {
        std::cout << "Error while parsing json file: " << errors << std::endl;
    }

    for (const auto &item : root)
    {
        int truckID = item["TruckID"].asInt();
        FleetData &data = FleetDataMap[truckID];

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

    int i = 1;

    for (ScalingTechnique ST : {ScalingTechnique::FIXEDAUTO, ScalingTechnique::FLEXIBLEAUTO, ScalingTechnique::FLEXIBLEAUTOEXT})
    {
        for (SecurityLevel SL : {SecurityLevel::HEStd_128_classic, SecurityLevel::HEStd_192_classic, SecurityLevel::HEStd_256_classic, SecurityLevel::HEStd_128_quantum, SecurityLevel::HEStd_192_quantum, SecurityLevel::HEStd_256_quantum})
        {
            for (int BS : {32})
            {
                for (auto ModSize : {std::make_pair(50, 49), std::make_pair(55, 50), std::make_pair(60, 29), std::make_pair(60, 45), std::make_pair(60, 50), std::make_pair(60, 59)}) // {std::make_pair(60, 22) removed as Error is too high}
                { //{std::make_pair(60, 22), std::make_pair(60, 59), std::make_pair(60, 40), std::make_pair(50, 22), std::make_pair(30, 22) }) {
                    std::cout << "##############################" << std::endl;
                    std::cout << "  Benchmark no " << i++ << "/270" << std::endl;
                    std::cout << "##############################" << std::endl;
                    Benchmark(Filename, Repetitions, ST, SL, BS, ModSize.first, ModSize.second);
                }
            }
        }
    }

    int milliseconds = TOC(t);
    int seconds = milliseconds / 1000;
    int hours = seconds / 3600;
    int remainingSeconds = seconds % 3600;
    int minutes = remainingSeconds / 60;

    std::cout << "\n Benchmarking completed in " << hours << " hours and " << minutes << " minutes." << std::endl;

    return 0;
}
