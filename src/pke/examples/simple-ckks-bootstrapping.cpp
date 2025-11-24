//==================================================================================
// BSD 2-Clause License
//
// Copyright (c) 2014-2022, NJIT, Duality Technologies Inc. and other contributors
//
// All rights reserved.
//
// Author TPOC: contact@openfhe.org
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//==================================================================================

/*

Example for CKKS bootstrapping with full packing

*/

#include "openfhe.h"

using namespace lbcrypto;

void SimpleBootstrapExample(usint dcrtBits, usint ringDim, usint numSlots, ScalingTechnique rescaleTech);
void SimpleBootstrapStCExample(usint dcrtBits, usint ringDim, usint numSlots, ScalingTechnique rescaleTech);

int main(int argc, char* argv[]) {
    uint32_t dcrtBits = 59;
    uint32_t ringDim  = 1 << 7;
    uint32_t numSlots = ringDim / 4;
    // SimpleBootstrapExample(dcrtBits, ringDim, numSlots, FIXEDMANUAL);
    // SimpleBootstrapStCExample(dcrtBits, ringDim, numSlots, FIXEDMANUAL);
    // SimpleBootstrapExample(dcrtBits, ringDim, numSlots, FIXEDAUTO);
    // SimpleBootstrapStCExample(dcrtBits, ringDim, numSlots, FIXEDAUTO);
    // SimpleBootstrapExample(dcrtBits, ringDim, numSlots, FLEXIBLEAUTO);
    // SimpleBootstrapStCExample(dcrtBits, ringDim, numSlots, FLEXIBLEAUTO);
    SimpleBootstrapExample(dcrtBits, ringDim, numSlots, FLEXIBLEAUTOEXT);
    SimpleBootstrapStCExample(dcrtBits, ringDim, numSlots, FLEXIBLEAUTOEXT);
}

void SimpleBootstrapExample(usint dcrtBits, usint ringDim, usint numSlots, ScalingTechnique rescaleTech) {
    std::cerr << "*****SimpleBootstrapExample with dcrtBits = " << dcrtBits << ", ringDim = " << ringDim
              << ", numSlots = " << numSlots << ", rescaleTech = " << rescaleTech << "*****" << std::endl;

    CCParams<CryptoContextCKKSRNS> parameters;
    // A. Specify main parameters
    /*  A1) Secret key distribution
    * The secret key distribution for CKKS should either be SPARSE_TERNARY or UNIFORM_TERNARY.
    * The SPARSE_TERNARY distribution was used in the original CKKS paper,
    * but in this example, we use UNIFORM_TERNARY because this is included in the homomorphic
    * encryption standard.
    */
    SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);

    /*  A2) Desired security level based on FHE standards.
    * In this example, we use the "NotSet" option, so the example can run more quickly with
    * a smaller ring dimension. Note that this should be used only in
    * non-production environments, or by experts who understand the security
    * implications of their choices. In production-like environments, we recommend using
    * HEStd_128_classic, HEStd_192_classic, or HEStd_256_classic for 128-bit, 192-bit,
    * or 256-bit security, respectively. If you choose one of these as your security level,
    * you do not need to set the ring dimension.
    */
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(ringDim);
    // parameters.SetBatchSize(numSlots);

    /*  A3) Scaling parameters.
    * By default, we set the modulus sizes and rescaling technique to the following values
    * to obtain a good precision and performance tradeoff. We recommend keeping the parameters
    * below unless you are an FHE expert.
    */
#if NATIVEINT == 128
    ScalingTechnique rescaleTech = FIXEDAUTO;
    uint32_t dcrtBits            = 78;
    uint32_t firstMod            = 89;
#else
    // ScalingTechnique rescaleTech = FIXEDMANUAL;
    // uint32_t dcrtBits            = dcrtBits;
    uint32_t firstMod = 60;
#endif

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);

    /*  A4) Multiplicative depth.
    * The goal of bootstrapping is to increase the number of available levels we have, or in other words,
    * to dynamically increase the multiplicative depth. However, the bootstrapping procedure itself
    * needs to consume a few levels to run. We compute the number of bootstrapping levels required
    * using GetBootstrapDepth, and add it to levelsAvailableAfterBootstrap to set our initial multiplicative
    * depth. We recommend using the input parameters below to get started.
    */
    std::vector<uint32_t> levelBudget = {2, 2};

    // Note that the actual number of levels avalailable after bootstrapping before next bootstrapping
    // will be levelsAvailableAfterBootstrap - 1 because an additional level
    // is used for scaling the ciphertext before next bootstrapping (in 64-bit CKKS bootstrapping)
    uint32_t levelsAvailableAfterBootstrap = 10;
    uint32_t depth = levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth(levelBudget, secretKeyDist);
    parameters.SetMultiplicativeDepth(depth);

    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);

    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
    cryptoContext->Enable(ADVANCEDSHE);
    cryptoContext->Enable(FHE);

    std::cout << "CKKS scheme is using ring dimension " << ringDim << " and depth " << depth << std::endl << std::endl;

    cryptoContext->EvalBootstrapSetup(levelBudget, {0, 0}, numSlots);

    auto keyPair = cryptoContext->KeyGen();
    cryptoContext->EvalMultKeyGen(keyPair.secretKey);
    cryptoContext->EvalBootstrapKeyGen(keyPair.secretKey, numSlots);

    std::vector<double> x = {-4.0, -3.0, -2.0, -1.0, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 4.0, 5.0};
    // std::vector<double> x = {-40, -30, -20, -10, -75, -5, -25, 0, 25, 5, 75, 10, 20, 30, 40, 50};

    if (x.size() < numSlots)
        x = Fill<double>(x, numSlots);

    // std::vector<double> x;
    // std::random_device rd;
    // std::mt19937 gen(rd());
    // std::uniform_real_distribution<> dis(0.0, 1.0);
    // for (size_t i = 0; i < numSlots; i++) {
    //     x.push_back(dis(gen));
    // }
    size_t encodedLength = x.size();

    // We start with a depleted ciphertext that has used up all of its levels.
    Plaintext ptxt = cryptoContext->MakeCKKSPackedPlaintext(x, 1, depth - 1, nullptr, numSlots);

    ptxt->SetLength(encodedLength);
    // std::cout << "Input: " << ptxt << "\n";

    std::cerr << "First numSlots elements of the input: [";
    std::copy_n(x.begin(), numSlots, std::ostream_iterator<double>(std::cerr, " "));
    std::cerr << "]" << std::endl;

    Ciphertext<DCRTPoly> ciph = cryptoContext->Encrypt(keyPair.publicKey, ptxt);

    std::cout << "Initial number of levels remaining: " << depth - ciph->GetLevel() << "\n\n";
    std::cerr << "Number of levels consumed by initial ciphertext: " << ciph->GetLevel()
              << ", and noise scale degree: " << ciph->GetNoiseScaleDeg() << std::endl;

    // ciph = cryptoContext->EvalMult(ciph, ciph);
    // std::cerr << "Number of levels consumed by ciphertext after multiplication: " << ciph->GetLevel()
    //           << ", and noise scale degree: " << ciph->GetNoiseScaleDeg() << std::endl;

    // auto start = std::chrono::high_resolution_clock::now();

    // Perform the bootstrapping operation. The goal is to increase the number of levels remaining
    // for HE computation.
    auto ciphertextAfter = cryptoContext->EvalBootstrap(ciph);

    // auto stop = std::chrono::high_resolution_clock::now();
    // std::cout << "Bootstrapping time: " << std::chrono::duration<double>(stop - start).count() << " s\n\n";

    std::cout << "Number of levels remaining after bootstrapping: "
              << depth - ciphertextAfter->GetLevel() - (ciphertextAfter->GetNoiseScaleDeg() - 1) << "\n\n";

    Plaintext result;
    cryptoContext->Decrypt(keyPair.secretKey, ciphertextAfter, &result);
    result->SetLength(encodedLength);
    // std::cout << "Output after bootstrapping: " << result << "\n";

    auto precision = result->GetRealPackedValue();

    std::cerr << "First numSlots elements of the output: [";
    std::copy_n(precision.begin(), numSlots, std::ostream_iterator<double>(std::cerr, " "));
    std::cerr << "]" << std::endl;

    std::transform(precision.begin(), precision.end(), x.begin(), precision.begin(), std::minus<double>());
    std::transform(precision.begin(), precision.end(), precision.begin(),
                   [&](const double& elem) { return std::floor(-log2(std::abs(elem))); });
    auto max_prec_it = std::max_element(precision.begin(), precision.end());
    auto min_prec_it = std::min_element(precision.begin(), precision.end());
    // std::cerr << "Decrypted precision: " << precision << "; max = " << *max_prec_it << ", min = " << *min_prec_it << std::endl << std::endl;
    std::cerr << "First numSlots elements of the decrypted precision: [";
    std::copy_n(precision.begin(), numSlots, std::ostream_iterator<double>(std::cerr, " "));
    std::cerr << "]; max = " << *max_prec_it << ", min = " << *min_prec_it << std::endl << std::endl;

    cryptoContext->ClearStaticMapsAndVectors();
}

void SimpleBootstrapStCExample(usint dcrtBits, usint ringDim, usint numSlots, ScalingTechnique rescaleTech) {
    std::cerr << "*****SimpleBootstrapStCExample with dcrtBits = " << dcrtBits << ", ringDim = " << ringDim
              << ", numSlots = " << numSlots << ", rescaleTech = " << rescaleTech << "*****" << std::endl;

    CCParams<CryptoContextCKKSRNS> parameters;
    // A. Specify main parameters
    /*  A1) Secret key distribution
    * The secret key distribution for CKKS should either be SPARSE_TERNARY or UNIFORM_TERNARY.
    * The SPARSE_TERNARY distribution was used in the original CKKS paper,
    * but in this example, we use UNIFORM_TERNARY because this is included in the homomorphic
    * encryption standard.
    */
    SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
    parameters.SetSecretKeyDist(secretKeyDist);

    /*  A2) Desired security level based on FHE standards.
    * In this example, we use the "NotSet" option, so the example can run more quickly with
    * a smaller ring dimension. Note that this should be used only in
    * non-production environments, or by experts who understand the security
    * implications of their choices. In production-like environments, we recommend using
    * HEStd_128_classic, HEStd_192_classic, or HEStd_256_classic for 128-bit, 192-bit,
    * or 256-bit security, respectively. If you choose one of these as your security level,
    * you do not need to set the ring dimension.
    */
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetRingDim(ringDim);
    // parameters.SetBatchSize(numSlots);

    /*  A3) Scaling parameters.
    * By default, we set the modulus sizes and rescaling technique to the following values
    * to obtain a good precision and performance tradeoff. We recommend keeping the parameters
    * below unless you are an FHE expert.
    */
#if NATIVEINT == 128
    ScalingTechnique rescaleTech = FIXEDAUTO;
    usint dcrtBits               = 78;
    usint firstMod               = 89;
#else
    // ScalingTechnique rescaleTech = FIXEDMANUAL;
    // usint dcrtBits               = dcrtBits;
    usint firstMod = 60;
#endif

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);

    /*  A4) Multiplicative depth.
    * The goal of bootstrapping is to increase the number of available levels we have, or in other words,
    * to dynamically increase the multiplicative depth. However, the bootstrapping procedure itself
    * needs to consume a few levels to run. We compute the number of bootstrapping levels required
    * using GetBootstrapDepth, and add it to levelsAvailableAfterBootstrap to set our initial multiplicative
    * depth. We recommend using the input parameters below to get started.
    */
    std::vector<uint32_t> levelBudget = {2, 2};

    // Note that the actual number of levels avalailable after bootstrapping before next bootstrapping
    // will be levelsAvailableAfterBootstrap - levelBudget[1] - 1 = 9 below because we perform SlotsToCoefficients
    // and a scaling the ciphertext before the modulus raise in the next bootstrapping (in 64-bit CKKS bootstrapping)
    uint32_t levelsAvailableAfterBootstrap = 10 + levelBudget[1];
    usint depth = levelsAvailableAfterBootstrap + FHECKKSRNS::GetBootstrapDepth({levelBudget[0], 0}, secretKeyDist);
    parameters.SetMultiplicativeDepth(depth);

    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);

    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
    cryptoContext->Enable(ADVANCEDSHE);
    cryptoContext->Enable(FHE);

    std::cout << "CKKS scheme is using ring dimension " << ringDim << " and depth " << depth << std::endl << std::endl;

    cryptoContext->EvalBootstrapSetup(levelBudget, {0, 0}, numSlots, 0, true, true);

    auto keyPair = cryptoContext->KeyGen();
    cryptoContext->EvalMultKeyGen(keyPair.secretKey);
    cryptoContext->EvalBootstrapKeyGen(keyPair.secretKey, numSlots);

    std::vector<double> x = {-4.0, -3.0, -2.0, -1.0, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 4.0, 5.0};
    // std::vector<double> x = {-40, -30, -20, -10, -75, -5, -25, 0, 25, 5, 75, 10, 20, 30, 40, 50};
    if (x.size() < numSlots)
        x = Fill<double>(x, numSlots);
    size_t encodedLength = x.size();

    // We start with a depleted ciphertext that has used up all of its levels, but which still
    // has sufficient levels to perform the SlotsToCoefficients operation.
    Plaintext ptxt = cryptoContext->MakeCKKSPackedPlaintext(x, 1, depth - 1 - levelBudget[1], nullptr, numSlots);
    // Plaintext ptxt = cryptoContext->MakeCKKSPackedPlaintext(
    //     x, 1, depth - 2 - levelBudget[1]); // Subtract and extra level to allow one multiplication before bootstrapping

    ptxt->SetLength(encodedLength);
    // std::cout << "Input: " << ptxt << std::endl;
    std::cerr << "First numSlots elements of the output: [";
    std::copy_n(x.begin(), numSlots, std::ostream_iterator<double>(std::cerr, " "));
    std::cerr << "]" << std::endl;

    Ciphertext<DCRTPoly> ciph = cryptoContext->Encrypt(keyPair.publicKey, ptxt);

    std::cout << "Initial number of levels remaining: " << depth - ciph->GetLevel() << " and total depth: " << depth
              << std::endl;
    std::cerr << "Number of levels consumed by initial ciphertext: " << ciph->GetLevel()
              << ", and noise scale degree: " << ciph->GetNoiseScaleDeg() << std::endl;

    // ciph = cryptoContext->EvalMult(ciph, ciph);
    // // cryptoContext->ModReduceInPlace(ciph);
    // std::cerr << "Number of levels consumed by ciphertext after multiplication: " << ciph->GetLevel()
    //           << ", and noise scale degree: " << ciph->GetNoiseScaleDeg() << std::endl;

    // Perform the bootstrapping operation. The goal is to increase the number of levels remaining
    // for HE computation.
    auto ciphertextAfter = cryptoContext->EvalBootstrap(ciph);

    std::cout << "Number of levels remaining after bootstrapping: "
              << depth - ciphertextAfter->GetLevel() - (ciphertextAfter->GetNoiseScaleDeg() - 1) << "\n\n";
    std::cerr << "Number of levels consumed by ciphertext after bootstrapping: " << ciphertextAfter->GetLevel()
              << ", and noise scale degree: " << ciphertextAfter->GetNoiseScaleDeg() << std::endl;

    Plaintext result;
    cryptoContext->Decrypt(keyPair.secretKey, ciphertextAfter, &result);
    result->SetLength(encodedLength);
    // std::cout << "Output after bootstrapping: " << result << "\n";

    auto precision = result->GetRealPackedValue();
    std::cerr << "First numSlots elements of the output: [";
    std::copy_n(precision.begin(), numSlots, std::ostream_iterator<double>(std::cerr, " "));
    std::cerr << "]" << std::endl;

    std::transform(precision.begin(), precision.end(), x.begin(), precision.begin(), std::minus<double>());
    std::transform(precision.begin(), precision.end(), precision.begin(),
                   [&](const double& elem) { return std::floor(-log2(std::abs(elem))); });
    auto max_prec_it = std::max_element(precision.begin(), precision.end());
    auto min_prec_it = std::min_element(precision.begin(), precision.end());
    // std::cerr << "Decrypted precision: " << precision << "; max = " << *max_prec_it << ", min = " << *min_prec_it << std::endl << std::endl;
    std::cerr << "First numSlots elements of the decrypted precision: [";
    std::copy_n(precision.begin(), numSlots, std::ostream_iterator<double>(std::cerr, " "));
    std::cerr << "]; max = " << *max_prec_it << ", min = " << *min_prec_it << std::endl << std::endl;

    cryptoContext->ClearStaticMapsAndVectors();
}
