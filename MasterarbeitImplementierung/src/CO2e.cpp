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
    double DurationCalculationCO2e;
    double DurationDecryption;
    int LevelForCO2e;

    std::string Comment = "";

    double SizeOfCiphertext;
    double ThroughputEncryption;
    double ThroughputDecryption;
    double ThroughputCO2eCalculation;

    int SizeCryptoContext;
    int SizeEvalMultKey;
    int SizeRotKey;
    int PublicKey;
    int SecretKey;

    std::string ActualCO2eResult;
    std::string ExpectedCO2eResult;
    double DurationCO2eCalculation;

    int PrecisionAfterCO2eCalculation;
};

class FHEBenchmark
{
private:
    CryptoContext<DCRTPoly> cryptoContext;
    KeyPair<DCRTPoly> keyPair;
    std::map<int, FleetDataCipher> CipherInput;
    std::vector<Ciphertext<DCRTPoly>> CO2eResult;
    uint32_t numSlots;
    usint depth;
    std::vector<double> RawCO2e;
    Results BenchmarkResults;

    std::vector<double> CalcPlaintextCO2e()
    {
        std::vector<double> co2eResults;
        for (const auto &[truckId, data] : FleetDataMap)
        {
            double co2e = 0.0;
            for (size_t i = 0; i < data.CargoWeight.size(); i++)
            {
                double cm1 = data.DistanceTravelled[i] * data.CargoWeight[i];
                double cm2 = data.EmissionFactor[i] * cm1;
                co2e += cm2;
            }
            // std::cout << "TruckID: " << truckId << ", CO2e: " << co2e << std::endl; // Debugging-Ausgabe
            co2eResults.push_back(co2e);
        }
        return co2eResults;
    }

    double CalculateAveragePrecisionBits(const std::vector<double> &result, const std::vector<double> &expectedResult)
    {
        if (result.size() != expectedResult.size())
            throw std::invalid_argument("Cannot compare vectors with different numbers of elements");

        double sumPrecisionBits = 0;
        size_t count = result.size();

        for (size_t i = 0; i < count; ++i)
        {
            // Berechne den relativen Fehler
            double relativeError = std::abs(result[i] - expectedResult[i]) / std::abs(expectedResult[i]);

            // Berechne die Präzision in Bits
            double precisionBits = -std::log2(relativeError);

            // Summe der Präzisionen in Bits
            sumPrecisionBits += precisionBits;
        }

        // Berechne den Mittelwert der Präzision in Bits
        double averagePrecisionBits = sumPrecisionBits / count;

        std::cout << "Average Precision (in Bits): " << averagePrecisionBits << std::endl;
        return averagePrecisionBits;
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
        if (SecLevel == SecurityLevel::HEStd_NotSet) // SecurityLevel::HEStd_NotSet
        {
            parameters.SetRingDim(8192);
        }

        parameters.SetBatchSize(batchSize);
        parameters.SetMultiplicativeDepth(2);
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
        cryptoContext->EvalRotateKeyGen(keyPair.secretKey, {1, 2, 3, -1, -2, -3});
        cryptoContext->EvalSumKeyGen(keyPair.secretKey, keyPair.publicKey);

        RawCO2e = CalcPlaintextCO2e();

        BenchmarkResults.batchSize = batchSize;
        BenchmarkResults.RescaleTech = RescaleTech;
        BenchmarkResults.SecLevel = SecLevel;
        BenchmarkResults.FirstModSize = FirstModsize;
        BenchmarkResults.ScaleModSize = ScaleModSize;
        BenchmarkResults.RingDimension = cryptoContext->GetRingDimension();
        BenchmarkResults.MultiplicativeDepth = 2;
        double log2q = log2(cryptoContext->GetCryptoParameters()->GetElementParams()->GetModulus().ConvertToDouble());
        BenchmarkResults.Log2q = std::isinf(log2q) ? 10000.0 : log2q;

        std::stringstream s1, s2, s3, s4, s5;
        size_t lengthInBytes;


        Serial::Serialize(cryptoContext, s1, SerType::BINARY);
        std::string str = s1.str();
        BenchmarkResults.SizeCryptoContext = str.size();

        cryptoContext->SerializeEvalMultKey(s2, SerType::BINARY);
        str = s2.str();
        lengthInBytes = str.size();
        BenchmarkResults.SizeEvalMultKey = lengthInBytes / 1024;

        cryptoContext->SerializeEvalAutomorphismKey(s3, SerType::BINARY);
        str = s3.str();
        lengthInBytes = str.size();
        BenchmarkResults.SizeRotKey = lengthInBytes / 1024;

        Serial::Serialize(keyPair.publicKey, s4, SerType::BINARY);
        str = s4.str();
        lengthInBytes = str.size();                                             // Korrekte Länge des Strings in Bytes
        BenchmarkResults.PublicKey = static_cast<double>(lengthInBytes) / 1024; // Größe in Kilobytes (KB)

        // Serialisierung und Größenberechnung des Secret Keys
        Serial::Serialize(keyPair.secretKey, s5, SerType::BINARY);
        str = s5.str();
        lengthInBytes = str.size();                                             // Korrekte Länge des Strings in Bytes
        BenchmarkResults.SecretKey = static_cast<double>(lengthInBytes) / 1024; // Größe i
    }

    void EncryptData()
    {
        TimeVar tStart, t;
        std::cout << "Start of encryption" << std::endl;
        std::vector<double> result;
        double durEncryption = 0.0;
        TIC(tStart);

        for (const auto &[truckId, data] : FleetDataMap)
        {
            FleetDataCipher &datacipher = CipherInput[truckId];

            // Encrypt CargoWeight
            Plaintext ptxt = cryptoContext->MakeCKKSPackedPlaintext(data.CargoWeight);
            ptxt->SetLength(numSlots);
            TIC(t);
            datacipher.CargoWeight = cryptoContext->Encrypt(keyPair.publicKey, ptxt);
            durEncryption += TOC_MS(t);

            // Encrypt CO2Emissions
            ptxt = cryptoContext->MakeCKKSPackedPlaintext(data.CO2Emissions);
            ptxt->SetLength(numSlots);
            TIC(t);
            datacipher.CO2Emissions = cryptoContext->Encrypt(keyPair.publicKey, ptxt);
            durEncryption += TOC_MS(t);

            // Encrypt DistanceTravelled
            ptxt = cryptoContext->MakeCKKSPackedPlaintext(data.DistanceTravelled);
            ptxt->SetLength(numSlots);
            TIC(t);
            datacipher.DistanceTravelled = cryptoContext->Encrypt(keyPair.publicKey, ptxt);
            durEncryption += TOC_MS(t);

            // Encrypt DriverID
            ptxt = cryptoContext->MakeCKKSPackedPlaintext(data.DriverID);
            ptxt->SetLength(numSlots);
            TIC(t);
            datacipher.DriverID = cryptoContext->Encrypt(keyPair.publicKey, ptxt);
            durEncryption += TOC_MS(t);

            // Encrypt EmissionFactor
            ptxt = cryptoContext->MakeCKKSPackedPlaintext(data.EmissionFactor);
            ptxt->SetLength(numSlots);
            TIC(t);
            datacipher.EmissionFactor = cryptoContext->Encrypt(keyPair.publicKey, ptxt);
            durEncryption += TOC_MS(t);

            // Encrypt FuelConsumption
            ptxt = cryptoContext->MakeCKKSPackedPlaintext(data.FuelConsumption);
            ptxt->SetLength(numSlots);
            TIC(t);
            datacipher.FuelConsumption = cryptoContext->Encrypt(keyPair.publicKey, ptxt);
            durEncryption += TOC_MS(t);

            // Encrypt LocLat
            ptxt = cryptoContext->MakeCKKSPackedPlaintext(data.LocLat);
            ptxt->SetLength(numSlots);
            TIC(t);
            datacipher.LocLat = cryptoContext->Encrypt(keyPair.publicKey, ptxt);
            durEncryption += TOC_MS(t);

            // Encrypt LocLong
            ptxt = cryptoContext->MakeCKKSPackedPlaintext(data.LocLong);
            ptxt->SetLength(numSlots);
            TIC(t);
            datacipher.LocLong = cryptoContext->Encrypt(keyPair.publicKey, ptxt);
            durEncryption += TOC_MS(t);

            // Encrypt Speed
            ptxt = cryptoContext->MakeCKKSPackedPlaintext(data.Speed);
            ptxt->SetLength(numSlots);
            TIC(t);
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

    void CalcCipherCO2e()
    {
        TimeVar t, tDecrypt;                                     // Timer variable to measure execution time
        TIC(t);                                                  // Start the timer
        std::cout << "\nStart of CO2e calculation" << std::endl; // Print a message indicating the start of the CO2e calculation

        CO2eResult.clear(); // Clear the previous results
        BenchmarkResults.LevelForCO2e = CipherInput.begin()->second.FuelConsumption->GetLevel();
        std::cout << "\tLevel before CO2e calculation: " << CipherInput.begin()->second.FuelConsumption->GetLevel() << std::endl;

        for (const auto &[truckId, data] : CipherInput)
        {
            auto cm1 = cryptoContext->EvalMult(data.DistanceTravelled, data.CargoWeight);
            auto cm2 = cryptoContext->EvalMult(data.EmissionFactor, cm1);
            cm2 = cryptoContext->EvalSum(cm2, numSlots); // Summing across all slots

            CO2eResult.push_back(cm2);
        }

        std::cout << "\tLevel after CO2e calculation: " << CO2eResult[0]->GetLevel() << std::endl;
        BenchmarkResults.LevelForCO2e = CO2eResult[0]->GetLevel();
        BenchmarkResults.DurationCalculationCO2e = round(TOC_MS(t) * 100) / 100.0;

        std::cout << "Total calculation duration: " << BenchmarkResults.DurationCalculationCO2e << " ms" << std::endl;

        try
        {
            TIC(tDecrypt); // Start the decryption timer
            // Decrypt the encrypted result
            std::vector<double> results;

            for (size_t i = 0; i < CO2eResult.size(); i++)
            {
                Plaintext resultCO2e;
                cryptoContext->Decrypt(keyPair.secretKey, CO2eResult[i], &resultCO2e);
                resultCO2e->SetLength(numSlots);
                results.push_back(resultCO2e->GetRealPackedValue()[0]); // All elements of plaintext are equal therefore its sufficient to read first element
            }

            // Record the precision of the CO2e result
            BenchmarkResults.PrecisionAfterCO2eCalculation = std::floor(CalculateAveragePrecisionBits(RawCO2e, results));
            BenchmarkResults.ExpectedCO2eResult = vectorToFormattedString(RawCO2e, RawCO2e.size());
            BenchmarkResults.ActualCO2eResult = vectorToFormattedString(results, results.size());

            // Print the results of the CO2e calculation
            std::cout << "\nResults of CO2e calculation" << std::endl;
            printVector("\tExpected result:\t", RawCO2e, RawCO2e.size());
            printVector("\tActual result:\t\t", results, results.size());
            std::cout << "\tPrecision: " << BenchmarkResults.PrecisionAfterCO2eCalculation << " bits" << std::endl;

            BenchmarkResults.DurationDecryption = round(TOC_MS(tDecrypt) * 100) / 100.0; // Record the decryption time
            std::cout << "Total decryption duration: " << BenchmarkResults.DurationDecryption << " ms" << std::endl;
        }
        catch (const std::exception &e)
        {
            // Handle decryption errors
            BenchmarkResults.Comment = "Error while decrypting CO2e results";
            std::cerr << e.what() << '\n';
        }
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
        double throughputCO2eCalculation = (BenchmarkResults.SizeOfCiphertext * length) / (static_cast<double>(BenchmarkResults.DurationCalculationCO2e) / 1000);

        BenchmarkResults.ThroughputEncryption = RoundToDecimals(throughputEncryption, 2);
        BenchmarkResults.ThroughputDecryption = RoundToDecimals(throughputDecryption, 2);
        BenchmarkResults.ThroughputCO2eCalculation = RoundToDecimals(throughputCO2eCalculation, 2);
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

        newResult["DurationCalculationCO2e [ms]"] = BenchmarkResults.DurationCalculationCO2e;
        newResult["PrecisionAfterCO2eCalculation [bit]"] = BenchmarkResults.PrecisionAfterCO2eCalculation;
        newResult["ActualCO2eResult"] = BenchmarkResults.ActualCO2eResult;
        newResult["ExpectedCO2eResult"] = BenchmarkResults.ExpectedCO2eResult;
        newResult["Duration Decryption [ms]"] = BenchmarkResults.DurationDecryption;
        newResult["Duration Encryption [ms]"] = BenchmarkResults.DurationEncryption;

        newResult["SizeCryptoContext [byte]"] = BenchmarkResults.SizeCryptoContext;
        newResult["SizeOfCiphertext [kb]"] = BenchmarkResults.SizeOfCiphertext;
        newResult["SizeEvalMultKey [kb]"] = BenchmarkResults.SizeEvalMultKey;
        newResult["SizeRotKey [kb]"] = BenchmarkResults.SizeRotKey;
        newResult["SizePublicKey [kb]"] = BenchmarkResults.PublicKey;
        newResult["SizeSecretKey[kb]"] = BenchmarkResults.SecretKey;

        newResult["ThroughputDecryption [kb/s]"] = RoundToDecimals(BenchmarkResults.ThroughputDecryption, 2);
        newResult["ThroughputEncryption [kb/s]"] = RoundToDecimals(BenchmarkResults.ThroughputEncryption, 2);
        newResult["ThroughputCO2eCalculation [kb/s]"] = RoundToDecimals(BenchmarkResults.ThroughputCO2eCalculation, 2);

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
        bm.CalcCipherCO2e();
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
        std::cerr << "Repetitions must be a positive integer" << std::endl;
        return 1;
    }

    Json::Value root;
    Json::CharReaderBuilder builder;
    std::ifstream file("FleetData.json");
    std::string errors;

    if (!Json::parseFromStream(builder, file, &root, &errors))
    {
        std::cout << "Error while parsing JSON file: " << errors << std::endl;
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

    // int i = 1;

    // for (ScalingTechnique ST : {ScalingTechnique::FIXEDAUTO, ScalingTechnique::FLEXIBLEAUTO, ScalingTechnique::FLEXIBLEAUTOEXT})
    // {
    //     for (SecurityLevel SL : {SecurityLevel::HEStd_128_classic, SecurityLevel::HEStd_192_classic, SecurityLevel::HEStd_256_classic, SecurityLevel::HEStd_128_quantum, SecurityLevel::HEStd_192_quantum, SecurityLevel::HEStd_256_quantum})
    //     {
    //         for (int BS : {32})
    //         {
    //             for (auto ModSize : {std::make_pair(55, 39)})//std::make_pair(50, 49), std::make_pair(55, 50), std::make_pair(60, 29), std::make_pair(60, 45), std::make_pair(60, 50), std::make_pair(60, 59)})
    //             {
    //                 std::cout << "##############################" << std::endl;
    //                 std::cout << "  Benchmark no " << i++ << std::endl;
    //                 std::cout << "##############################" << std::endl;
    //                 Benchmark(Filename, Repetitions, ST, SL, BS, ModSize.first, ModSize.second);
    //             }
    //         }
    //     }
    // }
    Benchmark(Filename, Repetitions, ScalingTechnique::FLEXIBLEAUTO, SecurityLevel::HEStd_192_quantum, 32, 60, 59);
    // Benchmark(Filename, Repetitions, ScalingTechnique::FLEXIBLEAUTO, SecurityLevel::HEStd_192_classic, 32, 60, 37);
    //  Benchmark(Filename, Repetitions, ScalingTechnique::FLEXIBLEAUTO, SecurityLevel::HEStd_256_quantum, 32, 60, 35);

    // Benchmark(Filename, Repetitions, ScalingTechnique::FLEXIBLEAUTO, SecurityLevel::HEStd_192_quantum, 32, 60, 30); //39

    int milliseconds = TOC(t);
    int seconds = milliseconds / 1000;
    int hours = seconds / 3600;
    int remainingSeconds = seconds % 3600;
    int minutes = remainingSeconds / 60;

    std::cout << "\nBenchmarking completed in " << hours << " hours and " << minutes << " minutes." << std::endl;

    return 0;
}
