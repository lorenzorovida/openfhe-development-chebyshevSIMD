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
CKKS implementation. See https://eprint.iacr.org/2020/1118 for details.
 */

#include "cryptocontext.h"
#include "scheme/ckksrns/ckksrns-cryptoparameters.h"
#include "scheme/ckksrns/ckksrns-advancedshe.h"
#include "scheme/ckksrns/ckksrns-utils.h"
#include "schemebase/base-scheme.h"

#include <complex>
#include <vector>

namespace lbcrypto {

Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalMultMany(const std::vector<Ciphertext<DCRTPoly>>& ciphertextVec,
                                                      const std::vector<EvalKey<DCRTPoly>>& evalKeys) const {
    const uint32_t inSize = ciphertextVec.size();

    if (inSize == 0)
        OPENFHE_THROW("Input ciphertext vector is empty.");

    if (inSize == 1)
        return ciphertextVec[0]->Clone();

    const uint32_t lim = inSize * 2 - 2;
    std::vector<Ciphertext<DCRTPoly>> ciphertextMultVec(inSize - 1);

    auto algo               = ciphertextVec[0]->GetCryptoContext()->GetScheme();
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersRNS>(ciphertextVec[0]->GetCryptoParameters());
    uint32_t levelsToDrop   = cryptoParams->GetCompositeDegree();

    uint32_t i = 0, j = 0;
    for (; i < (inSize - 1); i += 2) {
        ciphertextMultVec[j] = algo->EvalMultAndRelinearize(ciphertextVec[i], ciphertextVec[i + 1], evalKeys);
        algo->ModReduceInPlace(ciphertextMultVec[j++], levelsToDrop);
    }
    if (i < inSize) {
        ciphertextMultVec[j] =
            algo->EvalMultAndRelinearize(ciphertextVec[i], ciphertextMultVec[i + 1 - inSize], evalKeys);
        algo->ModReduceInPlace(ciphertextMultVec[j++], levelsToDrop);
        i += 2;
    }
    for (; i < lim; i += 2) {
        ciphertextMultVec[j] =
            algo->EvalMultAndRelinearize(ciphertextMultVec[i - inSize], ciphertextMultVec[i + 1 - inSize], evalKeys);
        algo->ModReduceInPlace(ciphertextMultVec[j++], levelsToDrop);
    }

    return ciphertextMultVec.back();
}

//------------------------------------------------------------------------------
// LINEAR WEIGHTED SUM
//------------------------------------------------------------------------------

template <typename VectorDataType>
Ciphertext<DCRTPoly> internalEvalLinearWSum(const std::vector<ReadOnlyCiphertext<DCRTPoly>>& ciphertexts,
                                            const std::vector<VectorDataType>& constants) {
    const uint32_t limit = ciphertexts.size();
    std::vector<Ciphertext<DCRTPoly>> cts(limit);
    for (uint32_t i = 0; i < limit; ++i)
        cts[i] = ciphertexts[i]->Clone();
    return internalEvalLinearWSumMutable(cts, constants);
}

template <typename VectorDataType>
static inline Ciphertext<DCRTPoly> internalEvalLinearWSumBatch(std::vector<ReadOnlyCiphertext<DCRTPoly>>& ciphertexts,
                                                          const std::vector<std::vector<VectorDataType>>& batchOfConstants) {
    std::vector<Ciphertext<DCRTPoly>> cts(ciphertexts.size());
    for (uint32_t i = 0; i < ciphertexts.size(); i++)
        cts[i] = ciphertexts[i]->Clone();
    return internalEvalLinearWSumMutableBatch(cts, batchOfConstants);
}

template <typename VectorDataType>
Ciphertext<DCRTPoly> internalEvalLinearWSumMutable(std::vector<Ciphertext<DCRTPoly>>& ciphertexts,
                                                   const std::vector<VectorDataType>& constants) {
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(ciphertexts[0]->GetCryptoParameters());

    auto cc = ciphertexts[0]->GetCryptoContext();

    const uint32_t limit = ciphertexts.size();

    if (cryptoParams->GetScalingTechnique() != FIXEDMANUAL) {
        // Check to see if input ciphertexts are of same level
        // and adjust if needed to the max level among them
        uint32_t maxLevel = ciphertexts[0]->GetLevel();
        uint32_t maxIdx   = 0;
        for (uint32_t i = 1; i < limit; ++i) {
            if ((ciphertexts[i]->GetLevel() > maxLevel) ||
                ((ciphertexts[i]->GetLevel() == maxLevel) && (ciphertexts[i]->GetNoiseScaleDeg() == 2))) {
                maxLevel = ciphertexts[i]->GetLevel();
                maxIdx   = i;
            }
        }

        auto algo = cc->GetScheme();
        for (uint32_t i = 0; i < maxIdx; ++i)
            algo->AdjustLevelsAndDepthInPlace(ciphertexts[i], ciphertexts[maxIdx]);
        for (uint32_t i = maxIdx + 1; i < limit; ++i)
            algo->AdjustLevelsAndDepthInPlace(ciphertexts[i], ciphertexts[maxIdx]);

        uint32_t compositeDegree = cryptoParams->GetCompositeDegree();
        if (ciphertexts[maxIdx]->GetNoiseScaleDeg() == 2) {
            for (uint32_t i = 0; i < limit; ++i)
                algo->ModReduceInternalInPlace(ciphertexts[i], compositeDegree);
        }
    }

    cc->EvalMultInPlace(ciphertexts[0], constants[0]);
    for (uint32_t i = 1; i < limit; ++i) {
        cc->EvalMultInPlace(ciphertexts[i], constants[i]);
        cc->EvalAddInPlaceNoCheck(ciphertexts[0], ciphertexts[i]);
    }
    cc->ModReduceInPlace(ciphertexts[0]);
    return ciphertexts[0];
}

template <typename VectorDataType>
Ciphertext<DCRTPoly> EvalPartialLinearWSum(const std::vector<Ciphertext<DCRTPoly>>& ciphertexts,
                                           const std::vector<VectorDataType>& constants, uint32_t limit = 0) {
    if (0 == limit)
        limit = ciphertexts.size();

    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(ciphertexts[0]->GetCryptoParameters());

    auto cc = ciphertexts[0]->GetCryptoContext();

    std::vector<Ciphertext<DCRTPoly>> cts(limit);
    if (cryptoParams->GetScalingTechnique() != FIXEDMANUAL) {
        cts[0] = ciphertexts[0]->Clone();
        // Check to see if input ciphertexts are of same level
        // and adjust if needed to the max level among them
        uint32_t maxLevel = cts[0]->GetLevel();
        uint32_t maxIdx   = 0;
        for (uint32_t i = 1; i < limit; ++i) {
            cts[i] = ciphertexts[i]->Clone();
            if ((cts[i]->GetLevel() > maxLevel) ||
                ((cts[i]->GetLevel() == maxLevel) && (cts[i]->GetNoiseScaleDeg() == 2))) {
                maxLevel = cts[i]->GetLevel();
                maxIdx   = i;
            }
        }

        auto algo = cc->GetScheme();
        auto& ctm = cts[maxIdx];
        for (uint32_t i = 0; i < maxIdx; ++i)
            algo->AdjustLevelsAndDepthInPlace(cts[i], ctm);
        for (uint32_t i = maxIdx + 1; i < limit; ++i)
            algo->AdjustLevelsAndDepthInPlace(cts[i], ctm);

        uint32_t compositeDegree = cryptoParams->GetCompositeDegree();
        if (ctm->GetNoiseScaleDeg() == 2) {
            for (uint32_t i = 0; i < limit; ++i)
                algo->ModReduceInternalInPlace(cts[i], compositeDegree);
        }
    }
    else {
        for (uint32_t i = 0; i < limit; ++i)
            cts[i] = ciphertexts[i]->Clone();
    }

    cc->EvalMultInPlace(cts[0], constants[1]);
    for (uint32_t i = 1; i < limit; ++i) {
        cc->EvalMultInPlace(cts[i], constants[i + 1]);
        cc->EvalAddInPlaceNoCheck(cts[0], cts[i]);
    }
    cc->ModReduceInPlace(cts[0]);
    return cts[0];
}

template <typename VectorDataType>
static inline Ciphertext<DCRTPoly> internalEvalLinearWSumMutableBatch(
    std::vector<Ciphertext<DCRTPoly>>& ciphertexts, const std::vector<std::vector<VectorDataType>>& setOfConstants) {
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(ciphertexts[0]->GetCryptoParameters());

    auto cc = ciphertexts[0]->GetCryptoContext();

    //std::cout << "Ciao dalla somma:)" << std::endl;
    //std::cout << ciphertexts.size() << std::endl;
    //std::cout << ciphertexts[0]->GetSlots() << std::endl;
    //std::cout << setOfConstants.size() << std::endl;
    //std::cout << setOfConstants[0].size() << std::endl;

    int repetitions = ciphertexts[0]->GetSlots() / setOfConstants.size();

    if (cryptoParams->GetScalingTechnique() != FIXEDMANUAL) {
        // Check to see if input ciphertexts are of same level
        // and adjust if needed to the max level among them
        uint32_t maxLevel = ciphertexts[0]->GetLevel();
        uint32_t maxIdx   = 0;
        for (uint32_t i = 1; i < ciphertexts.size(); ++i) {
            if ((ciphertexts[i]->GetLevel() > maxLevel) ||
                ((ciphertexts[i]->GetLevel() == maxLevel) && (ciphertexts[i]->GetNoiseScaleDeg() == 2))) {
                maxLevel = ciphertexts[i]->GetLevel();
                maxIdx   = i;
            }
        }

        auto algo = cc->GetScheme();
        for (uint32_t i = 0; i < maxIdx; ++i)
            algo->AdjustLevelsAndDepthInPlace(ciphertexts[i], ciphertexts[maxIdx]);
        for (uint32_t i = maxIdx + 1; i < ciphertexts.size(); ++i)
            algo->AdjustLevelsAndDepthInPlace(ciphertexts[i], ciphertexts[maxIdx]);

        uint32_t compositeDegree = cryptoParams->GetCompositeDegree();
        if (ciphertexts[maxIdx]->GetNoiseScaleDeg() == 2) {
            for (uint32_t i = 0; i < ciphertexts.size(); ++i) {
                algo->ModReduceInternalInPlace(ciphertexts[i], compositeDegree);
            }
        }
    }

    std::vector<Plaintext> constantsPtxt(setOfConstants[0].size());

    for (uint32_t i = 0; i < setOfConstants[0].size(); i++) {

        /*
         * Type check for CKKS encoding
         */

        std::vector<double> tmp;
        tmp.reserve(setOfConstants.size());

        for (uint32_t j = 0; j < setOfConstants.size(); j++) {
            if constexpr (std::is_same_v<VectorDataType, double>) {
                tmp.push_back(setOfConstants[j][i]);
            }
            else {
                OPENFHE_THROW(std::string("Not implemented for type: ") + typeid(VectorDataType).name());
            }
        }

        if (repetitions > 1) {
            std::vector<double> tmpRep;

            for (int i = 0; i < repetitions; i++) {
                tmpRep.insert(tmpRep.end(), tmp.begin(), tmp.end());
            }

            tmp = tmpRep;
        }

        //std::cout << i << "/" << setOfConstants[0].size() << ", new tmp: " << tmp.size() << ", " << tmp[0] << ", " << tmp[1] << ", " << tmp[2] << std::endl;

        constantsPtxt[i] = cc->MakeCKKSPackedPlaintext(tmp, 1, ciphertexts[0]->GetLevel(),
                                                       nullptr, ciphertexts[0]->GetSlots());

    }

    Ciphertext<DCRTPoly> weightedSum = cc->EvalMult(ciphertexts[0], constantsPtxt[0]);

    Ciphertext<DCRTPoly> tmp;
    for (uint32_t i = 1; i < ciphertexts.size(); i++) {
        tmp = cc->EvalMult(ciphertexts[i], constantsPtxt[i]);
        cc->EvalAddInPlace(weightedSum, tmp);
    }

    cc->ModReduceInPlace(weightedSum);

    //std::cout << "Me ne vado" << std::endl;

    return weightedSum;
}

Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalLinearWSum(std::vector<ReadOnlyCiphertext<DCRTPoly>>& ciphertexts,
                                                        const std::vector<int64_t>& constants) const {
    return internalEvalLinearWSum(ciphertexts, constants);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalLinearWSum(std::vector<ReadOnlyCiphertext<DCRTPoly>>& ciphertexts,
                                                        const std::vector<double>& constants) const {
    return internalEvalLinearWSum(ciphertexts, constants);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalLinearWSum(std::vector<ReadOnlyCiphertext<DCRTPoly>>& ciphertexts,
                                                        const std::vector<std::complex<double>>& constants) const {
    return internalEvalLinearWSum(ciphertexts, constants);
}

// For batched SIMD Chebyshev
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalLinearWSumBatch(
    std::vector<ReadOnlyCiphertext<DCRTPoly>>& ciphertexts,
    const std::vector<std::vector<int64_t>>& batchOfConstants) const {
    return internalEvalLinearWSumBatch(ciphertexts, batchOfConstants);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalLinearWSumBatch(
    std::vector<ReadOnlyCiphertext<DCRTPoly>>& ciphertexts,
    const std::vector<std::vector<double>>& batchOfConstants) const {
    return internalEvalLinearWSumBatch(ciphertexts, batchOfConstants);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalLinearWSumBatch(
    std::vector<ReadOnlyCiphertext<DCRTPoly>>& ciphertexts,
    const std::vector<std::vector<std::complex<double>>>& batchOfConstants) const {
    return internalEvalLinearWSumBatch(ciphertexts, batchOfConstants);
}


Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalLinearWSumMutable(std::vector<Ciphertext<DCRTPoly>>& ciphertexts,
                                                               const std::vector<int64_t>& constants) const {
    return internalEvalLinearWSumMutable(ciphertexts, constants);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalLinearWSumMutable(std::vector<Ciphertext<DCRTPoly>>& ciphertexts,
                                                               const std::vector<double>& constants) const {
    return internalEvalLinearWSumMutable(ciphertexts, constants);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalLinearWSumMutable(
    std::vector<Ciphertext<DCRTPoly>>& ciphertexts, const std::vector<std::complex<double>>& constants) const {
    return internalEvalLinearWSumMutable(ciphertexts, constants);
}

// For batched SIMD Chebyshev
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalLinearWSumMutableBatch(
    std::vector<Ciphertext<DCRTPoly>>& ciphertexts, const std::vector<std::vector<int64_t>>& batchOfConstants) const {
    return internalEvalLinearWSumMutableBatch(ciphertexts, batchOfConstants);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalLinearWSumMutableBatch(
    std::vector<Ciphertext<DCRTPoly>>& ciphertexts, const std::vector<std::vector<double>>& batchOfConstants) const {
    return internalEvalLinearWSumMutableBatch(ciphertexts, batchOfConstants);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalLinearWSumMutableBatch(
    std::vector<Ciphertext<DCRTPoly>>& ciphertexts,
    const std::vector<std::vector<std::complex<double>>>& batchOfConstants) const {
    return internalEvalLinearWSumMutableBatch(ciphertexts, batchOfConstants);
}

//------------------------------------------------------------------------------
// EVAL POLYNOMIAL
//------------------------------------------------------------------------------

template <typename VectorDataType>
std::shared_ptr<seriesPowers<DCRTPoly>> internalEvalPowersLinear(ConstCiphertext<DCRTPoly>& x,
                                                                 const std::vector<VectorDataType>& coefficients) {
    const uint32_t k = coefficients.size() - 1;
    std::vector<bool> indices(k);

    // find indices for powers of x that need to be computed
    for (uint32_t i = k; i > 0; --i) {
        if (0 == (i & (i - 1))) {  // if i is a power of 2
            indices[i - 1] = true;
        }
        else {  // non-power of 2
            if (IsNotEqualZero(coefficients[i])) {
                uint32_t rem = i;

                // while rem is not a power of 2
                // set indices required to compute rem to 1
                while (0 != (rem & (rem - 1))) {
                    indices[rem - 1] = true;
                    rem &= (uint64_t(1) << (GetMSB(rem) - 1)) - 1;
                }
            }
        }
    }

    std::vector<Ciphertext<DCRTPoly>> powers(k);
    powers[0] = x->Clone();
    auto cc   = x->GetCryptoContext();

    auto cryptoParams        = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(x->GetCryptoParameters());
    uint32_t compositeDegree = cryptoParams->GetCompositeDegree();

    // computes all powers up to k for x
    for (uint32_t i = 2; i <= k; ++i) {
        if (0 == (i & (i - 1))) {
            powers[i - 1] = cc->EvalSquare(powers[i / 2 - 1]);
            cc->ModReduceInPlace(powers[i - 1]);
        }
        else {
            if (indices[i - 1]) {
                uint64_t p    = (uint64_t(1) << (GetMSB(i) - 1)) - 1;
                uint64_t r    = (i & p) - 1;
                uint32_t diff = powers[p]->GetLevel() - powers[r]->GetLevel();
                cc->LevelReduceInPlace(powers[r], nullptr, diff / compositeDegree);

                powers[i - 1] = cc->EvalMult(powers[p], powers[r]);
                cc->ModReduceInPlace(powers[i - 1]);
            }
        }
    }

    // brings all powers of x to the same level
    for (uint32_t i = 1; i < k; ++i) {
        if (indices[i - 1]) {
            uint32_t diff = powers[k - 1]->GetLevel() - powers[i - 1]->GetLevel();
            cc->LevelReduceInPlace(powers[i - 1], nullptr, diff / compositeDegree);
        }
    }

    return std::make_shared<seriesPowers<DCRTPoly>>(std::move(powers));
}

std::shared_ptr<seriesPowers<DCRTPoly>> internalEvalPowersPS(ConstCiphertext<DCRTPoly>& x, uint32_t degree) {
    auto degs  = ComputeDegreesPS(degree);
    uint32_t k = degs[0];
    uint32_t m = degs[1];

    std::vector<Ciphertext<DCRTPoly>> powers(k);
    powers[0] = x->Clone();

    auto cc = x->GetCryptoContext();
    uint32_t compositeDegree =
        std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(x->GetCryptoParameters())->GetCompositeDegree();

    // computes all powers up to k for x
    uint32_t powerOf2 = 2;
    uint32_t rem      = 0;
    for (uint32_t i = 2; i <= k; ++i) {
        if (rem == 0) {
            powers[i - 1] = cc->EvalSquare(powers[(powerOf2 >> 1) - 1]);
        }
        else {
            uint32_t diff = powers[powerOf2 - 1]->GetLevel() - powers[rem - 1]->GetLevel();
            cc->LevelReduceInPlace(powers[rem - 1], nullptr, diff / compositeDegree);
            powers[i - 1] = cc->EvalMult(powers[powerOf2 - 1], powers[rem - 1]);
        }

        if (++rem == powerOf2) {
            powerOf2 <<= 1;
            rem = 0;
        }
        cc->ModReduceInPlace(powers[i - 1]);
    }

    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(powers[k - 1]->GetCryptoParameters());
    if (cryptoParams->GetScalingTechnique() == FIXEDMANUAL) {
        // brings all powers of x to the same level
        uint32_t levelk = powers[k - 1]->GetLevel();
        for (uint32_t i = 1; i < k; ++i)
            cc->LevelReduceInPlace(powers[i - 1], nullptr, levelk - powers[i - 1]->GetLevel());
    }
    else {
        for (uint32_t i = 1; i < k; ++i)
            cc->GetScheme()->AdjustLevelsAndDepthInPlace(powers[i - 1], powers[k - 1]);
    }

    // computes powers of form k*2^i for x and the product of the powers in power2, that yield x^{k(2*m - 1)}
    std::vector<Ciphertext<DCRTPoly>> powers2(m);
    powers2[0] = powers.back();

    auto power2km1 = powers.back();

    for (uint32_t i = 1; i < m; ++i) {
        powers2[i] = cc->EvalSquare(powers2[i - 1]);
        cc->ModReduceInPlace(powers2[i]);
        power2km1 = cc->EvalMult(powers2[i], power2km1);
        cc->ModReduceInPlace(power2km1);
    }

    return std::make_shared<seriesPowers<DCRTPoly>>(std::move(powers), std::move(powers2), std::move(power2km1), k, m);
}

std::shared_ptr<seriesPowers<DCRTPoly>> AdvancedSHECKKSRNS::EvalPowers(ConstCiphertext<DCRTPoly>& x,
                                                                       const std::vector<int64_t>& coefficients) const {
    uint32_t d = Degree(coefficients);
    return (d < 5) ? internalEvalPowersLinear(x, coefficients) : internalEvalPowersPS(x, d);
}
std::shared_ptr<seriesPowers<DCRTPoly>> AdvancedSHECKKSRNS::EvalPowers(ConstCiphertext<DCRTPoly>& x,
                                                                       const std::vector<double>& coefficients) const {
    uint32_t d = Degree(coefficients);
    return (d < 5) ? internalEvalPowersLinear(x, coefficients) : internalEvalPowersPS(x, d);
}
std::shared_ptr<seriesPowers<DCRTPoly>> AdvancedSHECKKSRNS::EvalPowers(
    ConstCiphertext<DCRTPoly>& x, const std::vector<std::complex<double>>& coefficients) const {
    uint32_t d = Degree(coefficients);
    return (d < 5) ? internalEvalPowersLinear(x, coefficients) : internalEvalPowersPS(x, d);
}

template <typename VectorDataType>
Ciphertext<DCRTPoly> internalEvalPolyLinearWithPrecomp(std::vector<Ciphertext<DCRTPoly>>& powers,
                                                       const std::vector<VectorDataType>& coefficients) {
    const uint32_t k = coefficients.size() - 1;
    if (k <= 1)
        OPENFHE_THROW("The coefficients vector should contain at least 2 elements");

    if (!IsNotEqualZero(coefficients[k]))
        OPENFHE_THROW("EvalPolyLinear: The highest-order coefficient cannot be set to 0.");

    auto cc = powers[0]->GetCryptoContext();

    // perform scalar multiplication for the highest-order term
    auto result = cc->EvalMult(powers[k - 1], coefficients[k]);

    // perform scalar multiplication for all other terms and sum them up
    for (uint32_t i = 0; i < k - 1; ++i) {
        if (IsNotEqualZero(coefficients[i + 1])) {
            cc->EvalMultInPlace(powers[i], coefficients[i + 1]);
            cc->EvalAddInPlace(result, powers[i]);
        }
    }

    // Do rescaling after scalar multiplication
    cc->ModReduceInPlace(result);

    // adds the free term (at x^0)
    cc->EvalAddInPlace(result, coefficients[0]);

    return result;
}

template <typename VectorDataType>
Ciphertext<DCRTPoly> InnerEvalPolyPS(ConstCiphertext<DCRTPoly>& x, const std::vector<VectorDataType>& coefficients,
                                     uint32_t k, uint32_t m, const std::vector<Ciphertext<DCRTPoly>>& powers,
                                     const std::vector<Ciphertext<DCRTPoly>>& powers2) {
    // Compute k*2^m because we use it often
    uint32_t k2m2k = k * (1 << (m - 1)) - k;

    // Divide coefficients by x^{k*2^{m-1}}
    std::vector<VectorDataType> xkm(k2m2k + k + 1);
    xkm.back() = 1;
    auto divqr = LongDivisionPoly(coefficients, xkm);

    // Subtract x^{k(2^{m-1} - 1)} from r
    auto& r2 = divqr->r;
    if (auto n = Degree(r2); static_cast<int32_t>(k2m2k - n) <= 0) {
        r2.resize(n + 1);
        r2[k2m2k] -= 1;
    }
    else {
        r2.resize(k2m2k + 1);
        r2.back() = -1;
    }

    auto divcs = LongDivisionPoly(r2, divqr->q);
    auto cc    = x->GetCryptoContext();

    Ciphertext<DCRTPoly> cu, qu, su;

#pragma omp task shared(qu)
    {
        // Evaluate q and s2 at u.
        // If their degrees are larger than k, then recursively apply the Paterson-Stockmeyer algorithm.

        if (Degree(divqr->q) > k) {
            qu = InnerEvalPolyPS(x, divqr->q, k, m - 1, powers, powers2);
        }
        else {
            qu = cc->EvalAdd(powers[k - 1], divqr->q.front());
            divqr->q.resize(k);
            if (uint32_t n = Degree(divqr->q); n > 0)
                cc->EvalAddInPlace(qu, EvalPartialLinearWSum(powers, divqr->q, n));
        }
    }

#pragma omp task shared(su)
    {
        // Add x^{k(2^{m-1} - 1)} to s
        auto& s2 = divcs->r;
        s2.resize(k2m2k + 1);
        s2.back() = 1;

        if (Degree(s2) > k) {
            su = InnerEvalPolyPS(x, s2, k, m - 1, powers, powers2);
        }
        else {
            su = cc->EvalAdd(powers[k - 1], s2.front());
            s2.resize(k);
            if (uint32_t n = Degree(s2); n > 0)
                cc->EvalAddInPlace(su, EvalPartialLinearWSum(powers, s2, n));
        }
    }

    if (uint32_t n = Degree(divcs->q); n == 0) {
        cu = cc->EvalAdd(powers2[m - 1], divcs->q.front());
    }
    else if (n == 1) {
        if (IsNotEqualOne(divcs->q[1])) {
            cu = cc->EvalMult(powers.front(), divcs->q[1]);
            cc->ModReduceInPlace(cu);
            cc->EvalAddInPlace(cu, powers2[m - 1]);
        }
        else {
            cu = cc->EvalAdd(powers2[m - 1], powers.front());
        }
        cc->EvalAddInPlace(cu, divcs->q.front());
    }
    else {
        cu = cc->EvalAdd(powers2[m - 1], EvalPartialLinearWSum(powers, divcs->q, n));
        cc->EvalAddInPlace(cu, divcs->q.front());
    }

#pragma omp taskwait

    auto result = cc->EvalMult(cu, qu);
    cc->ModReduceInPlace(result);
    cc->EvalAddInPlace(result, su);
    return result;
}

template <typename VectorDataType>
Ciphertext<DCRTPoly> internalEvalPolyPSWithPrecomp(const std::shared_ptr<seriesPowers<DCRTPoly>>& ctxtPowers,
                                                   const std::vector<VectorDataType>& coefficients) {
    auto& powers    = ctxtPowers->powersRe;
    auto& powers2   = ctxtPowers->powers2Re;
    auto& power2km1 = ctxtPowers->power2km1Re;
    auto k          = ctxtPowers->k;
    auto m          = ctxtPowers->m;

    // Compute k*2^{m-1}-k because we use it a lot
    uint32_t k2m2k = k * (1 << (m - 1)) - k;

    // Add T^{k(2^m - 1)}(y) to the polynomial that has to be evaluated
    auto f2 = coefficients;
    f2.resize(Degree(f2) + 1);
    f2.resize(2 * k2m2k + k + 1);
    f2.back() = 1;

    Ciphertext<DCRTPoly> result;
#pragma omp parallel num_threads(OpenFHEParallelControls.GetThreadLimit(6 * m + 2))
    {
#pragma omp single
        result =
            powers[0]->GetCryptoContext()->EvalSub(InnerEvalPolyPS(powers[0], f2, k, m, powers, powers2), power2km1);
    }
    return result;
}

Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalPoly(ConstCiphertext<DCRTPoly>& x,
                                                  const std::vector<int64_t>& coeffs) const {
    return (Degree(coeffs) < 5) ? EvalPolyLinear(x, coeffs) : EvalPolyPS(x, coeffs);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalPoly(ConstCiphertext<DCRTPoly>& x,
                                                  const std::vector<double>& coeffs) const {
    return (Degree(coeffs) < 5) ? EvalPolyLinear(x, coeffs) : EvalPolyPS(x, coeffs);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalPoly(ConstCiphertext<DCRTPoly>& x,
                                                  const std::vector<std::complex<double>>& coeffs) const {
    return (Degree(coeffs) < 5) ? EvalPolyLinear(x, coeffs) : EvalPolyPS(x, coeffs);
}

Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalPolyWithPrecomp(std::shared_ptr<seriesPowers<DCRTPoly>> ctxtPowers,
                                                             const std::vector<int64_t>& coeffs) const {
    return (Degree(coeffs) < 5) ? internalEvalPolyLinearWithPrecomp(ctxtPowers->powersRe, coeffs) :
                                  internalEvalPolyPSWithPrecomp(ctxtPowers, coeffs);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalPolyWithPrecomp(std::shared_ptr<seriesPowers<DCRTPoly>> ctxtPowers,
                                                             const std::vector<double>& coeffs) const {
    return (Degree(coeffs) < 5) ? internalEvalPolyLinearWithPrecomp(ctxtPowers->powersRe, coeffs) :
                                  internalEvalPolyPSWithPrecomp(ctxtPowers, coeffs);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalPolyWithPrecomp(std::shared_ptr<seriesPowers<DCRTPoly>> ctxtPowers,
                                                             const std::vector<std::complex<double>>& coeffs) const {
    return (Degree(coeffs) < 5) ? internalEvalPolyLinearWithPrecomp(ctxtPowers->powersRe, coeffs) :
                                  internalEvalPolyPSWithPrecomp(ctxtPowers, coeffs);
}

Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalPolyLinear(ConstCiphertext<DCRTPoly>& x,
                                                        const std::vector<int64_t>& coeffs) const {
    return internalEvalPolyLinearWithPrecomp(internalEvalPowersLinear(x, coeffs)->powersRe, coeffs);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalPolyLinear(ConstCiphertext<DCRTPoly>& x,
                                                        const std::vector<double>& coeffs) const {
    return internalEvalPolyLinearWithPrecomp(internalEvalPowersLinear(x, coeffs)->powersRe, coeffs);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalPolyLinear(ConstCiphertext<DCRTPoly>& x,
                                                        const std::vector<std::complex<double>>& coeffs) const {
    return internalEvalPolyLinearWithPrecomp(internalEvalPowersLinear(x, coeffs)->powersRe, coeffs);
}

Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalPolyPS(ConstCiphertext<DCRTPoly>& x,
                                                    const std::vector<int64_t>& coeffs) const {
    return internalEvalPolyPSWithPrecomp(internalEvalPowersPS(x, Degree(coeffs)), coeffs);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalPolyPS(ConstCiphertext<DCRTPoly>& x,
                                                    const std::vector<double>& coeffs) const {
    return internalEvalPolyPSWithPrecomp(internalEvalPowersPS(x, Degree(coeffs)), coeffs);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalPolyPS(ConstCiphertext<DCRTPoly>& x,
                                                    const std::vector<std::complex<double>>& coeffs) const {
    return internalEvalPolyPSWithPrecomp(internalEvalPowersPS(x, Degree(coeffs)), coeffs);
}

//------------------------------------------------------------------------------
// EVAL CHEBYSHEV SERIES
//------------------------------------------------------------------------------

template <typename VectorDataType>
std::shared_ptr<seriesPowers<DCRTPoly>> internalEvalChebyPolysLinear(ConstCiphertext<DCRTPoly>& x,
                                                                     const std::vector<VectorDataType>& coefficients,
                                                                     double a, double b) {
    const uint32_t k = coefficients.size() - 1;
    std::vector<Ciphertext<DCRTPoly>> T(k);

    auto cc = x->GetCryptoContext();

    // computes linear transformation y = -1 + 2 (x-a)/(b-a)
    // consumes one level when a <> -1 && b <> 1
    if (!IsNotEqualNegOne(a) && !IsNotEqualOne(b)) {
        T[0] = x->Clone();
    }
    else {
        // linear transformation is needed
        double alpha = 2 / (b - a);
        double beta  = a * alpha;

        T[0] = cc->EvalMult(x, alpha);
        cc->ModReduceInPlace(T[0]);
        cc->EvalAddInPlace(T[0], -1.0 - beta);
    }

    // Computes Chebyshev polynomials up to degree k
    // for y: T_1(y) = y, T_2(y), ... , T_k(y)
    // uses binary tree multiplication
    for (uint32_t i = 2; i <= k; ++i) {
        if (i & 0x1) {  // if i is odd
            // compute T_{2i+1}(y) = 2*T_i(y)*T_{i+1}(y) - y
            T[i - 1] = cc->EvalMult(T[i / 2 - 1], T[i / 2]);
            cc->EvalAddInPlaceNoCheck(T[i - 1], T[i - 1]);
            cc->ModReduceInPlace(T[i - 1]);
            cc->EvalSubInPlace(T[i - 1], T[0]);
        }
        else {
            // compute T_{2i}(y) = 2*T_i(y)^2 - 1
            T[i - 1] = cc->EvalSquare(T[i / 2 - 1]);
            cc->EvalAddInPlaceNoCheck(T[i - 1], T[i - 1]);
            cc->ModReduceInPlace(T[i - 1]);
            cc->EvalAddInPlace(T[i - 1], -1.0);
        }
    }

    uint32_t compositeDegree =
        std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(x->GetCryptoParameters())->GetCompositeDegree();
    for (uint32_t i = 1; i < k; ++i)
        cc->LevelReduceInPlace(T[i - 1], nullptr, (T[k - 1]->GetLevel() - T[i - 1]->GetLevel()) / compositeDegree);

    return std::make_shared<seriesPowers<DCRTPoly>>(std::move(T));
}

template <typename VectorDataType>
Ciphertext<DCRTPoly> internalEvalChebyshevSeriesLinearWithPrecomp(std::vector<Ciphertext<DCRTPoly>>& T,
                                                                  const std::vector<VectorDataType>& coefficients) {
    const uint32_t k = coefficients.size() - 2;

    // perform scalar multiplication for the highest-order term
    auto cc     = T[0]->GetCryptoContext();
    auto result = cc->EvalMult(T[k], coefficients[k + 1]);

    // perform scalar multiplication for all other terms and sum them up
    for (uint32_t i = 0; i < k; ++i) {
        if (IsNotEqualZero(coefficients[i + 1])) {
            cc->EvalMultInPlace(T[i], coefficients[i + 1]);
            cc->EvalAddInPlace(result, T[i]);
        }
    }

    // Do rescaling after scalar multiplication
    cc->ModReduceInPlace(result);

    // adds the free term (at x^0)
    cc->EvalAddInPlace(result, coefficients[0] / 2.0);

    return result;
}

template <typename VectorDataType>
Ciphertext<DCRTPoly> InnerEvalChebyshevPS(ConstCiphertext<DCRTPoly>& x, const std::vector<VectorDataType>& coefficients,
                                          uint32_t k, uint32_t m, const std::vector<Ciphertext<DCRTPoly>>& T,
                                          const std::vector<Ciphertext<DCRTPoly>>& T2) {
    // Compute k*2^{m-1}-k because we use it a lot
    uint32_t k2m2k = k * (1 << (m - 1)) - k;

    // Divide coefficients by T^{k*2^{m-1}}
    std::vector<VectorDataType> Tkm(k2m2k + k + 1);
    Tkm.back() = 1;
    auto divqr = LongDivisionChebyshev(coefficients, Tkm);

    // Subtract x^{k(2^{m-1} - 1)} from r
    auto& r2 = divqr->r;
    if (uint32_t n = Degree(r2); static_cast<int32_t>(k2m2k - n) <= 0) {
        r2.resize(n + 1);
        r2[k2m2k] -= 1;
    }
    else {
        r2.resize(k2m2k + 1);
        r2.back() = -1;
    }

    auto divcs = LongDivisionChebyshev(r2, divqr->q);
    auto cc    = x->GetCryptoContext();

    Ciphertext<DCRTPoly> cu, qu, su;

    {
        // Evaluate q and s2 at u.
        // If their degrees are larger than k, then recursively apply the Paterson-Stockmeyer algorithm.
        if (Degree(divqr->q) > k) {
            qu = InnerEvalChebyshevPS(x, divqr->q, k, m - 1, T, T2);
        }
        else {
            // dq = k from construction
            // perform scalar multiplication for all other terms and sum them up if there are non-zero coefficients

            // the highest order coefficient will always be a power of two up to 2^{m-1} because q is "monic" but the Chebyshev rule adds a factor of 2
            // we don't need to increase the depth by multiplying the highest order coefficient, but instead checking and summing, since we work with m <= 4.
            qu                   = T[k - 1]->Clone();
            const uint32_t limit = std::log2(ToReal(divqr->q.back()));
            for (uint32_t i = 0; i < limit; ++i)
                cc->EvalAddInPlaceNoCheck(qu, qu);

            // adds the free term (at x^0)
            cc->EvalAddInPlace(qu, divqr->q.front() / 2.0);
            // The number of levels of qu is the same as the number of levels of T[k-1] + 1.
            // Will only get here when m = 2, so the number of levels of qu and T2[m-1] will be the same.

            divqr->q.resize(k);
            if (uint32_t n = Degree(divqr->q); n > 0)
                cc->EvalAddInPlace(qu, EvalPartialLinearWSum(T, divqr->q, n));
        }
    }

    {
        // Add x^{k(2^{m-1} - 1)} to s
        auto& s2 = divcs->r;
        s2.resize(k2m2k + 1);
        s2.back() = 1;

        if (Degree(s2) > k) {
            su = InnerEvalChebyshevPS(x, s2, k, m - 1, T, T2);
        }
        else {
            // the highest order coefficient will always be 1 because s2 is monic.
            su = T[k - 1]->Clone();

            // ds = k from construction
            // perform scalar multiplication for all other terms and sum them up if there are non-zero coefficients
            s2.resize(k);
            if (uint32_t n = Degree(s2); n > 0)
                cc->EvalAddInPlace(su, EvalPartialLinearWSum(T, s2, n));

            // adds the free term (at x^0)
            cc->EvalAddInPlace(su, s2.front() / 2.0);

            // The number of levels of su is the same as the number of levels of T[k-1] or T[k-1] + 1. Need to reduce it to T2[m-1] + 1.
            cc->LevelReduceInPlace(su, nullptr);
        }
    }

    if (uint32_t n = Degree(divcs->q); n >= 1) {
        if (n == 1) {
            if (IsNotEqualOne(divcs->q[1])) {
                cu = cc->EvalMult(T.front(), divcs->q[1]);
                cc->ModReduceInPlace(cu);
            }
            else {
                cu = T.front()->Clone();
            }
        }
        else {
            cu = EvalPartialLinearWSum(T, divcs->q, n);
        }

        // adds the free term (at x^0)
        cc->EvalAddInPlace(cu, divcs->q.front() / 2.0);

        // Need to reduce levels up to the level of T2[m-1].
        uint32_t cd =
            std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(x->GetCryptoParameters())->GetCompositeDegree();
        cc->LevelReduceInPlace(cu, nullptr, (T2[m - 1]->GetLevel() - cu->GetLevel()) / cd);
    }

    cu = cu ? cc->EvalAdd(T2[m - 1], cu) : cc->EvalAdd(T2[m - 1], divcs->q.front() / 2.0);

    auto result = cc->EvalMult(cu, qu);
    cc->ModReduceInPlace(result);
    cc->EvalAddInPlace(result, su);
    return result;
}

template <typename VectorDataType>
static Ciphertext<DCRTPoly> InnerEvalChebyshevPSBatch(ConstCiphertext<DCRTPoly>& x,
                                                 const std::vector<std::vector<VectorDataType>>& batchOfCoefficients,
                                                 uint32_t k, uint32_t m, std::vector<Ciphertext<DCRTPoly>>& T,
                                                 std::vector<Ciphertext<DCRTPoly>>& T2) {

    std::vector<std::vector<double>> batchOfCoefficientsDouble;
    batchOfCoefficientsDouble.reserve(batchOfCoefficients.size());

    // This is required as throughout the function we will create CKKS plaintexts,
    // but the encoding function requires doubles (and not complex values, accepted by
    // VectorDataType)
    if constexpr (std::is_same_v<VectorDataType, double>) {
        for (const auto& innerVec : batchOfCoefficients) {
            std::vector<double> converted;
            converted.reserve(innerVec.size());

            for (const auto& value : innerVec) {
                converted.push_back(static_cast<double>(value));
            }
            batchOfCoefficientsDouble.push_back(std::move(converted));
        }
    } else {
        OPENFHE_THROW(std::string("Not implemented for type: ") + typeid(VectorDataType).name());
    }

    std::vector<std::shared_ptr<longDiv<double>>> divcsVec;
    std::vector<std::shared_ptr<longDiv<double>>> divqrVec;
    std::vector<std::vector<double>> s2Vec;

    auto cc = x->GetCryptoContext();
    uint32_t compositeDegree =
        std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(x->GetCryptoParameters())->GetCompositeDegree();

    for (size_t i = 0; i < batchOfCoefficientsDouble.size(); i++) {
        // Compute k*2^{m-1}-k because we use it a lot
        uint32_t k2m2k = k * (1 << (m - 1)) - k;

        // Divide coefficients by T^{k*2^{m-1}}
        std::vector<VectorDataType> Tkm(static_cast<int32_t>(k2m2k + k) + 1);
        Tkm.back() = 1;
        auto divqr = LongDivisionChebyshev(batchOfCoefficientsDouble[i], Tkm);

        // Subtract x^{k(2^{m-1} - 1)} from r
        auto r2 = divqr->r;
        if (static_cast<int32_t>(k2m2k - Degree(divqr->r)) <= 0) {
            r2[static_cast<int32_t>(k2m2k)] -= 1;
            r2.resize(Degree(r2) + 1);
        }
        else {
            r2.resize(static_cast<int32_t>(k2m2k + 1));
            r2.back() = -1;
        }

        // Divide r2 by q
        auto divcs = LongDivisionChebyshev(r2, divqr->q);

        // Add x^{k(2^{m-1} - 1)} to s
        auto s2 = divcs->r;
        s2.resize(static_cast<int32_t>(k2m2k + 1), 0.0);
        s2.back() = 1;

        divqrVec.push_back(divqr);
        divcsVec.push_back(divcs);
        s2Vec.push_back(s2);
    }

    // Evaluate c at u
    Ciphertext<DCRTPoly> cu;
    uint32_t dc = Degree(divcsVec[0]->q);
    bool flag_c = false;
    if (dc >= 1) {
        if (dc == 1) {
            //if (IsNotEqualOne(divcs->q[1])) {
            // NOTE: We can not optimize anymore since it is not a constant, so we always perform the product

            std::vector<double> coeffs;
            for (size_t i = 0; i < divcsVec.size(); i++) coeffs.push_back(divcsVec[i]->q[1]);

            cu = cc->EvalMult(T.front(), cc->MakeCKKSPackedPlaintext(coeffs, 1, T.front()->GetLevel(),
                                                                     nullptr, T.front()->GetSlots()));
            cc->ModReduceInPlace(cu);
        }
        else {
            std::vector<std::vector<double>> batchOfWeights;
            std::vector<Ciphertext<DCRTPoly>> ctxs(dc);

            for (uint32_t i = 0; i < dc; i++) {
                ctxs[i] = T[i];
            }

            for (size_t j = 0; j < divcsVec.size(); j++) {
                std::vector<double> weights(dc);
                for (uint32_t i = 0; i < dc; i++) {
                    weights[i] = divcsVec[j]->q[i + 1];
                }

                batchOfWeights.push_back(weights);
            }

            cu = internalEvalLinearWSumMutableBatch(ctxs, batchOfWeights);

        }

        // adds the free term (at x^0)
        std::vector<double> freeTerms;
        for (size_t i = 0; i < divcsVec.size(); i++) {
            freeTerms.push_back(divcsVec[i]->q.front() / 2.0);
        }

        Plaintext freeTermsPtxt = cc->MakeCKKSPackedPlaintext(freeTerms, 1, cu->GetLevel(),
                                                              nullptr, cu->GetSlots());
        cu = cc->EvalAdd(cu, freeTermsPtxt);
        // Need to reduce levels up to the level of T2[m-1].
        uint32_t levelDiff = T2[m - 1]->GetLevel() - cu->GetLevel();
        cc->LevelReduceInPlace(cu, nullptr, levelDiff / compositeDegree);

        flag_c = true;
    }

    // Evaluate q and s2 at u. If their degrees are larger than k, then recursively apply the Paterson-Stockmeyer algorithm.
    Ciphertext<DCRTPoly> qu;

    if (Degree(divqrVec[0]->q) > k) {
        std::vector<std::vector<double>> coeffs;
        for (size_t i = 0; i < divqrVec.size(); i++) {
            coeffs.push_back(divqrVec[i]->q);
        }

        qu = InnerEvalChebyshevPSBatch(x, coeffs, k, m - 1, T, T2);
    }
    else {
        // dq = k from construction
        // perform scalar multiplication for all other terms and sum them up if there are non-zero coefficients

        auto qcopy = divqrVec[0]->q;

        qcopy.resize(k);


        if (Degree(qcopy) > 0) {
            std::vector<Ciphertext<DCRTPoly>> ctxs(Degree(qcopy));
            std::vector<std::vector<double>> weights(divqrVec.size());

            for (size_t j = 0; j < divqrVec.size(); j++) {
                std::vector<double> weightsForLevel;
                for (uint32_t i = 0; i < Degree(qcopy); ++i) {
                    ctxs[i] = T[i];
                    weightsForLevel.push_back(divqrVec[j]->q[i + 1]);
                }
                weights[j] = weightsForLevel;
            }

            qu = cc->EvalLinearWSumMutableBatch(ctxs, weights);

            // the highest order coefficient will always be a power of two up to 2^{m-1} because q is "monic" but the Chebyshev rule adds a factor of 2
            // we don't need to increase the depth by multiplying the highest order coefficient, but instead checking and summing, since we work with m <= 4.
            Ciphertext<DCRTPoly> sum = T[k - 1]->Clone();
            uint32_t limit           = log2(ToReal(divqrVec[0]->q.back()));
            for (uint32_t i = 0; i < limit; ++i) {
                sum = cc->EvalAdd(sum, sum);
            }
            cc->EvalAddInPlace(qu, sum);
        }
        else {
            Ciphertext<DCRTPoly> sum = T[k - 1]->Clone();
            uint32_t limit           = log2(ToReal(divqrVec[0]->q.back()));
            for (uint32_t i = 0; i < limit; ++i) {
                sum = cc->EvalAdd(sum, sum);
            }
            qu = sum;
        }

        // adds the free term (at x^0)
        std::vector<double> freeTerms;
        for (size_t i = 0; i < divqrVec.size(); i++) {
            freeTerms.push_back(divqrVec[i]->q.front() / 2.0);
        }

        Plaintext ptxt = cc->MakeCKKSPackedPlaintext(freeTerms, 1, qu->GetLevel(),
                                                     nullptr, qu->GetSlots());
        qu = cc->EvalAdd(qu, ptxt);
        // The number of levels of qu is the same as the number of levels of T[k-1] or T[k-1] + 1.
        // No need to reduce it to T2[m-1] because it only reaches here when m = 2.
    }

    Ciphertext<DCRTPoly> su;

    if (Degree(s2Vec[0]) > k) {
        su = InnerEvalChebyshevPSBatch(x, s2Vec, k, m - 1, T, T2);
    }
    else {
        // ds = k from construction
        // perform scalar multiplication for all other terms and sum them up if there are non-zero coefficients
        auto scopy = s2Vec[0];
        scopy.resize(k);
        if (Degree(scopy) > 0) {
            std::vector<Ciphertext<DCRTPoly>> ctxs(Degree(scopy));
            std::vector<std::vector<double>> weights(s2Vec.size());

            for (uint32_t j = 0; j < s2Vec.size(); j++) {
                std::vector<double> weightsForLevel;
                for (uint32_t i = 0; i < Degree(scopy); ++i) {
                    ctxs[i] = T[i];
                    weightsForLevel.push_back(s2Vec[j][i + 1]);
                }
                weights[j] = weightsForLevel;
            }

            su = cc->EvalLinearWSumMutableBatch(ctxs, weights);

            // the highest order coefficient will always be 1 because s2 is monic.
            cc->EvalAddInPlace(su, T[k - 1]);
        }
        else {
            su = T[k - 1]->Clone();
        }

        // adds the free term (at x^0)
        std::vector<double> freeTerms;
        for (uint32_t i = 0; i < s2Vec.size(); i++) {
            freeTerms.push_back(s2Vec[i].front() / 2.0);
        }

        Plaintext ptxt = cc->MakeCKKSPackedPlaintext(freeTerms, 1, su->GetLevel(), nullptr, su->GetSlots());
        cc->EvalAddInPlace(su, ptxt);

        // The number of levels of su is the same as the number of levels of T[k-1] or T[k-1] + 1. Need to reduce it to T2[m-1] + 1.
        // su = cc->LevelReduce(su, nullptr, su->GetElements()[0].GetNumOfElements() - Lm + 1) ;
        cc->LevelReduceInPlace(su, nullptr);
    }

    Ciphertext<DCRTPoly> result;

    if (flag_c) {
        result = cc->EvalAdd(T2[m - 1], cu);
    }
    else {
        std::vector<double> terms;
        for (size_t i = 0; i < divcsVec.size(); i++) {
            terms.push_back(divcsVec[i]->q.front() / 2.0);
        }

        Plaintext ptxt = cc->MakeCKKSPackedPlaintext(terms, 1, T2[m - 1]->GetLevel(), nullptr, T2[m - 1]->GetSlots());
        result = cc->EvalAdd(T2[m - 1], ptxt);
    }

    result = cc->EvalMult(result, qu);
    cc->ModReduceInPlace(result);

    cc->EvalAddInPlace(result, su);

    return result;
}

template <typename VectorDataType>
static Ciphertext<DCRTPoly> InnerEvalChebyshevPSBatchRepeated(ConstCiphertext<DCRTPoly>& x,
                                                      const std::vector<std::vector<VectorDataType>>& batchOfCoefficients,
                                                      uint32_t k, uint32_t m, std::vector<Ciphertext<DCRTPoly>>& T,
                                                      std::vector<Ciphertext<DCRTPoly>>& T2) {

    std::vector<std::vector<double>> batchOfCoefficientsDouble;
    batchOfCoefficientsDouble.reserve(batchOfCoefficients.size());


    // This is required as throughout the function we will create CKKS plaintexts,
    // but the encoding function requires doubles (and not complex values, accepted by
    // VectorDataType)
    if constexpr (std::is_same_v<VectorDataType, double>) {
        for (const auto& innerVec : batchOfCoefficients) {
            std::vector<double> converted;
            converted.reserve(innerVec.size());

            for (const auto& value : innerVec) {
                converted.push_back(static_cast<double>(value));
            }
            batchOfCoefficientsDouble.push_back(std::move(converted));
        }
    } else {
        OPENFHE_THROW(std::string("Not implemented for type: ") + typeid(VectorDataType).name());
    }

    std::vector<std::shared_ptr<longDiv<double>>> divcsVec;
    std::vector<std::shared_ptr<longDiv<double>>> divqrVec;
    std::vector<std::vector<double>> s2Vec;

    int repetitions = T.front()->GetSlots() / batchOfCoefficientsDouble.size();

    auto cc = x->GetCryptoContext();
    uint32_t compositeDegree =
        std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(x->GetCryptoParameters())->GetCompositeDegree();

    for (size_t i = 0; i < batchOfCoefficientsDouble.size(); i++) {
        // Compute k*2^{m-1}-k because we use it a lot
        uint32_t k2m2k = k * (1 << (m - 1)) - k;

        // Divide coefficients by T^{k*2^{m-1}}
        std::vector<VectorDataType> Tkm(static_cast<int32_t>(k2m2k + k) + 1);
        Tkm.back() = 1;
        auto divqr = LongDivisionChebyshev(batchOfCoefficientsDouble[i], Tkm);

        // Subtract x^{k(2^{m-1} - 1)} from r
        auto r2 = divqr->r;
        if (static_cast<int32_t>(k2m2k - Degree(divqr->r)) <= 0) {
            r2[static_cast<int32_t>(k2m2k)] -= 1;
            r2.resize(Degree(r2) + 1);
        }
        else {
            r2.resize(static_cast<int32_t>(k2m2k + 1));
            r2.back() = -1;
        }

        // Divide r2 by q
        auto divcs = LongDivisionChebyshev(r2, divqr->q);

        // Add x^{k(2^{m-1} - 1)} to s
        auto s2 = divcs->r;
        s2.resize(static_cast<int32_t>(k2m2k + 1), 0.0);
        s2.back() = 1;

        divqrVec.push_back(divqr);
        divcsVec.push_back(divcs);
        s2Vec.push_back(s2);
    }

    // Evaluate c at u
    Ciphertext<DCRTPoly> cu;
    uint32_t dc = Degree(divcsVec[0]->q);
    bool flag_c = false;
    if (dc >= 1) {
        if (dc == 1) {
            //if (IsNotEqualOne(divcs->q[1])) {
            // NOTE: We can not optimize anymore since it is not a constant, so we always perform the product

            std::vector<double> coeffs;
            for (size_t i = 0; i < divcsVec.size(); i++) coeffs.push_back(divcsVec[i]->q[1]);

            std::vector<double> coeffsRep;
            for (int i = 0; i < repetitions; ++i) {
                coeffsRep.insert(coeffsRep.end(), coeffs.begin(), coeffs.end());
            }

            cu = cc->EvalMult(T.front(), cc->MakeCKKSPackedPlaintext(coeffsRep, 1, T.front()->GetLevel(),
                                                                     nullptr, T.front()->GetSlots()));
            cc->ModReduceInPlace(cu);
        }
        else {
            std::vector<std::vector<double>> batchOfWeights;
            std::vector<Ciphertext<DCRTPoly>> ctxs(dc);

            for (uint32_t i = 0; i < dc; i++) {
                ctxs[i] = T[i];
            }

            for (size_t j = 0; j < divcsVec.size(); j++) {
                std::vector<double> weights(dc);
                for (uint32_t i = 0; i < dc; i++) {
                    weights[i] = divcsVec[j]->q[i + 1];
                }

                batchOfWeights.push_back(weights);
            }

            cu = internalEvalLinearWSumMutableBatch(ctxs, batchOfWeights);

        }

        // adds the free term (at x^0)
        std::vector<double> freeTerms;
        for (size_t i = 0; i < divcsVec.size(); i++) {
            freeTerms.push_back(divcsVec[i]->q.front() / 2.0);
        }

        std::vector<double> freeTermsRep;
        for (int i = 0; i < repetitions; ++i) {
            freeTermsRep.insert(freeTermsRep.end(), freeTerms.begin(), freeTerms.end());
        }

        Plaintext freeTermsPtxt = cc->MakeCKKSPackedPlaintext(freeTermsRep, 1, cu->GetLevel(),
                                                              nullptr, cu->GetSlots());
        cu = cc->EvalAdd(cu, freeTermsPtxt);
        // Need to reduce levels up to the level of T2[m-1].
        uint32_t levelDiff = T2[m - 1]->GetLevel() - cu->GetLevel();
        cc->LevelReduceInPlace(cu, nullptr, levelDiff / compositeDegree);

        flag_c = true;
    }

    // Evaluate q and s2 at u. If their degrees are larger than k, then recursively apply the Paterson-Stockmeyer algorithm.
    Ciphertext<DCRTPoly> qu;

    if (Degree(divqrVec[0]->q) > k) {
        std::vector<std::vector<double>> coeffs;
        for (size_t i = 0; i < divqrVec.size(); i++) {
            coeffs.push_back(divqrVec[i]->q);
        }

        qu = InnerEvalChebyshevPSBatchRepeated(x, coeffs, k, m - 1, T, T2);
    }
    else {
        // dq = k from construction
        // perform scalar multiplication for all other terms and sum them up if there are non-zero coefficients

        auto qcopy = divqrVec[0]->q;

        qcopy.resize(k);


        if (Degree(qcopy) > 0) {
            std::vector<Ciphertext<DCRTPoly>> ctxs(Degree(qcopy));
            std::vector<std::vector<double>> weights(divqrVec.size());

            for (size_t j = 0; j < divqrVec.size(); j++) {
                std::vector<double> weightsForLevel;
                for (uint32_t i = 0; i < Degree(qcopy); ++i) {
                    ctxs[i] = T[i];
                    weightsForLevel.push_back(divqrVec[j]->q[i + 1]);
                }
                weights[j] = weightsForLevel;
            }

            qu = cc->EvalLinearWSumMutableBatch(ctxs, weights);

            // the highest order coefficient will always be a power of two up to 2^{m-1} because q is "monic" but the Chebyshev rule adds a factor of 2
            // we don't need to increase the depth by multiplying the highest order coefficient, but instead checking and summing, since we work with m <= 4.
            Ciphertext<DCRTPoly> sum = T[k - 1]->Clone();
            uint32_t limit           = log2(ToReal(divqrVec[0]->q.back()));
            for (uint32_t i = 0; i < limit; ++i) {
                sum = cc->EvalAdd(sum, sum);
            }
            cc->EvalAddInPlace(qu, sum);
        }
        else {
            Ciphertext<DCRTPoly> sum = T[k - 1]->Clone();
            uint32_t limit           = log2(ToReal(divqrVec[0]->q.back()));
            for (uint32_t i = 0; i < limit; ++i) {
                sum = cc->EvalAdd(sum, sum);
            }
            qu = sum;
        }

        // adds the free term (at x^0)
        std::vector<double> freeTerms;
        for (size_t i = 0; i < divqrVec.size(); i++) {
            freeTerms.push_back(divqrVec[i]->q.front() / 2.0);
        }

        std::vector<double> freeTermsRep;
        for (int i = 0; i < repetitions; ++i) {
            freeTermsRep.insert(freeTermsRep.end(), freeTerms.begin(), freeTerms.end());
        }

        Plaintext ptxt = cc->MakeCKKSPackedPlaintext(freeTermsRep, 1, qu->GetLevel(),
                                                     nullptr, qu->GetSlots());
        qu = cc->EvalAdd(qu, ptxt);
        // The number of levels of qu is the same as the number of levels of T[k-1] or T[k-1] + 1.
        // No need to reduce it to T2[m-1] because it only reaches here when m = 2.
    }

    Ciphertext<DCRTPoly> su;

    if (Degree(s2Vec[0]) > k) {
        su = InnerEvalChebyshevPSBatchRepeated(x, s2Vec, k, m - 1, T, T2);
    }
    else {
        // ds = k from construction
        // perform scalar multiplication for all other terms and sum them up if there are non-zero coefficients
        auto scopy = s2Vec[0];
        scopy.resize(k);
        if (Degree(scopy) > 0) {
            std::vector<Ciphertext<DCRTPoly>> ctxs(Degree(scopy));
            std::vector<std::vector<double>> weights(s2Vec.size());

            for (uint32_t j = 0; j < s2Vec.size(); j++) {
                std::vector<double> weightsForLevel;
                for (uint32_t i = 0; i < Degree(scopy); ++i) {
                    ctxs[i] = T[i];
                    weightsForLevel.push_back(s2Vec[j][i + 1]);
                }
                weights[j] = weightsForLevel;
            }

            su = cc->EvalLinearWSumMutableBatch(ctxs, weights);

            // the highest order coefficient will always be 1 because s2 is monic.
            cc->EvalAddInPlace(su, T[k - 1]);
        }
        else {
            su = T[k - 1]->Clone();
        }

        // adds the free term (at x^0)
        std::vector<double> freeTerms;
        for (uint32_t i = 0; i < s2Vec.size(); i++) {
            freeTerms.push_back(s2Vec[i].front() / 2.0);
        }

        std::vector<double> freeTermsRep;
        for (int i = 0; i < repetitions; ++i) {
            freeTermsRep.insert(freeTermsRep.end(), freeTerms.begin(), freeTerms.end());
        }

        Plaintext ptxt = cc->MakeCKKSPackedPlaintext(freeTermsRep, 1, su->GetLevel(), nullptr, su->GetSlots());
        cc->EvalAddInPlace(su, ptxt);

        // The number of levels of su is the same as the number of levels of T[k-1] or T[k-1] + 1. Need to reduce it to T2[m-1] + 1.
        // su = cc->LevelReduce(su, nullptr, su->GetElements()[0].GetNumOfElements() - Lm + 1) ;
        cc->LevelReduceInPlace(su, nullptr);
    }

    Ciphertext<DCRTPoly> result;

    if (flag_c) {
        result = cc->EvalAdd(T2[m - 1], cu);
    }
    else {
        std::vector<double> terms;
        for (size_t i = 0; i < divcsVec.size(); i++) {
            terms.push_back(divcsVec[i]->q.front() / 2.0);
        }

        std::vector<double> termsRep;
        for (int i = 0; i < repetitions; ++i) {
            termsRep.insert(termsRep.end(), terms.begin(), terms.end());
        }

        Plaintext ptxt = cc->MakeCKKSPackedPlaintext(termsRep, 1, T2[m - 1]->GetLevel(), nullptr, T2[m - 1]->GetSlots());
        result = cc->EvalAdd(T2[m - 1], ptxt);
    }

    result = cc->EvalMult(result, qu);
    cc->ModReduceInPlace(result);

    cc->EvalAddInPlace(result, su);

    return result;
}

static std::shared_ptr<seriesPowers<DCRTPoly>> internalEvalChebyPolysPS(ConstCiphertext<DCRTPoly>& x, uint32_t degree,
                                                                 double a, double b) {
    auto degs  = ComputeDegreesPS(degree);
    uint32_t k = degs[0];
    uint32_t m = degs[1];

    // computes linear transformation y = -1 + 2 (x-a)/(b-a)
    // consumes one level when a <> -1 && b <> 1
    auto cc = x->GetCryptoContext();
    std::vector<Ciphertext<DCRTPoly>> T(k);
    if (!IsNotEqualNegOne(a) && !IsNotEqualOne(b)) {
        // no linear transformation is needed if a = -1, b = 1
        // T_1(y) = y
        T[0] = x->Clone();
    }
    else {
        // linear transformation is needed
        double alpha = 2 / (b - a);
        double beta  = a * alpha;

        T[0] = cc->EvalMult(x, alpha);
        cc->ModReduceInPlace(T[0]);
        cc->EvalAddInPlace(T[0], -1.0 - beta);
    }

    // Computes Chebyshev polynomials up to degree k
    // for y: T_1(y) = y, T_2(y), ... , T_k(y)
    // uses binary tree multiplication
    for (uint32_t i = 2; i <= k; ++i) {
        if (i & 0x1) {  // if i is odd
            // compute T_{2i+1}(y) = 2*T_i(y)*T_{i+1}(y) - y
            T[i - 1] = cc->EvalMult(T[i / 2 - 1], T[i / 2]);
            cc->EvalAddInPlaceNoCheck(T[i - 1], T[i - 1]);
            cc->ModReduceInPlace(T[i - 1]);
            cc->EvalSubInPlace(T[i - 1], T[0]);
        }
        else {
            // compute T_{2i}(y) = 2*T_i(y)^2 - 1
            T[i - 1] = cc->EvalSquare(T[i / 2 - 1]);
            cc->EvalAddInPlaceNoCheck(T[i - 1], T[i - 1]);
            cc->ModReduceInPlace(T[i - 1]);
            cc->EvalAddInPlace(T[i - 1], -1.0);
        }
    }

    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(T[k - 1]->GetCryptoParameters());
    if (cryptoParams->GetScalingTechnique() == FIXEDMANUAL) {
        // brings all powers of x to the same level
        for (uint32_t i = 1; i < k; ++i)
            cc->LevelReduceInPlace(T[i - 1], nullptr, T[k - 1]->GetLevel() - T[i - 1]->GetLevel());
    }
    else {
        for (uint32_t i = 1; i < k; ++i)
            cc->GetScheme()->AdjustLevelsAndDepthInPlace(T[i - 1], T[k - 1]);
    }

    std::vector<Ciphertext<DCRTPoly>> T2(m);
    // T2[0] is used as a placeholder
    T2[0] = T.back();

    // computes T_{k(2*m - 1)}(y)
    auto T2km1 = T.back();

    for (uint32_t i = 1; i < m; ++i) {
        // Compute the Chebyshev polynomials T_k(y), T_{2k}(y), T_{4k}(y), ... , T_{2^{m-1}k}(y)
        T2[i] = cc->EvalSquare(T2[i - 1]);
        cc->EvalAddInPlaceNoCheck(T2[i], T2[i]);
        cc->ModReduceInPlace(T2[i]);
        cc->EvalAddInPlace(T2[i], -1.0);

        // compute T_{k(2*m - 1)} = 2*T_{k(2^{m-1}-1)}(y)*T_{k*2^{m-1}}(y) - T_k(y)
        T2km1 = cc->EvalMult(T2km1, T2[i]);
        cc->EvalAddInPlaceNoCheck(T2km1, T2km1);
        cc->ModReduceInPlace(T2km1);
        cc->EvalSubInPlace(T2km1, T2[0]);
    }

    return std::make_shared<seriesPowers<DCRTPoly>>(std::move(T), std::move(T2), std::move(T2km1), k, m);
}

template <typename VectorDataType>
Ciphertext<DCRTPoly> internalEvalChebyshevSeriesPSWithPrecomp(const std::shared_ptr<seriesPowers<DCRTPoly>>& ctxtPolys,
                                                              const std::vector<VectorDataType>& coefficients) {
    auto& T     = ctxtPolys->powersRe;
    auto& T2    = ctxtPolys->powers2Re;
    auto& T2km1 = ctxtPolys->power2km1Re;
    auto k      = ctxtPolys->k;
    auto m      = ctxtPolys->m;

    // Compute k*2^{m-1}-k because we use it a lot
    uint32_t k2m2k = k * (1 << (m - 1)) - k;

    // Add T^{k(2^m - 1)}(y) to the polynomial that has to be evaluated
    auto f2 = coefficients;
    f2.resize(Degree(f2) + 1);
    f2.resize(2 * k2m2k + k + 1);
    f2.back() = 1;

    return T[0]->GetCryptoContext()->EvalSub(InnerEvalChebyshevPS(T[0], f2, k, m, T, T2), T2km1);
}

template <typename VectorDataType>
static inline Ciphertext<DCRTPoly> internalEvalChebyshevSeriesPSBatchWithPrecomp(
    std::shared_ptr<seriesPowers<DCRTPoly>> ctxtPolys,
    const std::vector<std::vector<VectorDataType>>& batchOfCoefficients) {

    std::vector<std::vector<double>> batchOfCoefficientsDouble;
    batchOfCoefficientsDouble.reserve(batchOfCoefficients.size());

    // This is required as throughout the function we will create CKKS plaintexts,
    // but the encoding function requires doubles (and not complex values, accepted by
    // VectorDataType)
    if constexpr (std::is_same_v<VectorDataType, double>) {
        for (const auto& innerVec : batchOfCoefficients) {
            std::vector<double> converted;
            converted.reserve(innerVec.size());

            for (const auto& value : innerVec) {
                converted.push_back(static_cast<double>(value));
            }
            batchOfCoefficientsDouble.push_back(std::move(converted));
        }
    } else {
        OPENFHE_THROW(std::string("Not implemented for type: ") + typeid(VectorDataType).name());
    }

    std::vector<std::shared_ptr<longDiv<double>>> divcsVec;
    std::vector<std::shared_ptr<longDiv<double>>> divqrVec;
    std::vector<std::vector<double>> s2Vec;

    /*
     * Check on the batch of coefficients: they must all be of the same size
     */
    std::size_t size = batchOfCoefficients[0].size();
    for (const auto& inner : batchOfCoefficients) {
        if (inner.size() != size) {
            OPENFHE_THROW("The batch must contain coefficients of the same size!");
        }
    }

    auto T     = ctxtPolys->powersRe;
    auto T2    = ctxtPolys->powers2Re;
    auto T2km1 = ctxtPolys->power2km1Re;
    auto k     = ctxtPolys->k;
    auto m     = ctxtPolys->m;

    // Compute k*2^{m-1}-k because we use it a lot
    uint32_t k2m2k = k * (1 << (m - 1)) - k;

    for (size_t b = 0; b < batchOfCoefficientsDouble.size(); b++) {
        auto f2 = batchOfCoefficientsDouble[b];
        auto n  = Degree(f2);
        f2.resize(n + 1);

        // Add T^{k(2^m - 1)}(y) to the polynomial that has to be evaluated
        f2.resize(2 * k2m2k + k + 1, 0.0);
        f2.back() = 1;

        // Divide f2 by T^{k*2^{m-1}}
        std::vector<double> Tkm(k2m2k + k + 1);
        Tkm.back() = 1;

        auto divqr = LongDivisionChebyshev(f2, Tkm);

        // Subtract x^{k(2^{m-1} - 1)} from r
        auto r2 = divqr->r;
        if (static_cast<int32_t>(k2m2k - Degree(r2)) <= 0) {
            r2[static_cast<int32_t>(k2m2k)] -= 1;
            r2.resize(Degree(r2) + 1);
        }
        else {
            r2.resize(static_cast<int32_t>(k2m2k + 1));
            r2.back() = -1;
        }

        // Divide r2 by q
        auto divcs = LongDivisionChebyshev(r2, divqr->q);

        // Add x^{k(2^{m-1} - 1)} to s
        auto s2 = divcs->r;
        s2.resize(k2m2k + 1);
        s2.back() = 1;

        divqrVec.push_back(divqr);
        divcsVec.push_back(divcs);
        s2Vec.push_back(s2);
    }

    auto cc = T[0]->GetCryptoContext();

    // Evaluate c at u
    Ciphertext<DCRTPoly> cu;

    // We use the degree of the first one, they should all be the same anyways
    uint32_t dc = Degree(divcsVec[0]->q);
    bool flag_c = false;
    if (dc >= 1) {
        if (dc == 1) {
            // if (IsNotEqualOne(divcs->q[1])) {
            // NOTE: We can not optimize anymore since it is not a constant, so we always perform the product
            std::vector<double> coeffs;
            for (size_t i = 0; i < divcsVec.size(); i++) coeffs.push_back(divcsVec[i]->q[1]);

            cu = cc->EvalMult(T.front(), cc->MakeCKKSPackedPlaintext(coeffs, 1, T.front()->GetLevel(),
                                                                     nullptr, T.front()->GetSlots()));
            cc->ModReduceInPlace(cu);
        }
        else {
            std::vector<std::vector<double>> batchOfWeights;
            std::vector<Ciphertext<DCRTPoly>> ctxs(dc);

            for (size_t i = 0; i < dc; i++) {
                ctxs[i] = T[i];
            }

            for (size_t j = 0; j < divcsVec.size(); j++) {
                std::vector<double> weights(dc);
                for (size_t i = 0; i < dc; i++) {
                    weights[i] = divcsVec[j]->q[i + 1];
                }

                batchOfWeights.push_back(weights);
            }

            cu = internalEvalLinearWSumMutableBatch(ctxs, batchOfWeights);
        }

        // adds the free term (at x^0)
        std::vector<double> freeTerms;
        for (size_t i = 0; i < divcsVec.size(); i++) {
            freeTerms.push_back(divcsVec[i]->q.front() / 2.0);
        }

        Plaintext freeTermsPtxt = cc->MakeCKKSPackedPlaintext(freeTerms, 1, cu->GetLevel(),
                                                                       nullptr, cu->GetSlots());
        cu = cc->EvalAdd(cu, freeTermsPtxt);

        // TODO : Andrey why not T2[m-1]->GetLevel() instead?
        // Need to reduce levels to the level of T2[m-1].
        //    uint32_t levelDiff = y->GetLevel() - cu->GetLevel() + ceil(log2(k)) + m - 1;
        //    cc->LevelReduceInPlace(cu, nullptr, levelDiff);

        flag_c = true;
    }

    // Evaluate q and s2 at u. If their degrees are larger than k, then recursively apply the Paterson-Stockmeyer algorithm.
    Ciphertext<DCRTPoly> qu;

    // Again, degrees should all be the same for the different divqrVec
    if (Degree(divqrVec[0]->q) > k) {
        std::vector<std::vector<double>> coeffs;
        for (size_t i = 0; i < divqrVec.size(); i++) {
            coeffs.push_back(divqrVec[i]->q);
        }
        qu = InnerEvalChebyshevPSBatch(T[0], coeffs, k, m - 1, T, T2);
    }
    else {
        // dq = k from construction
        // perform scalar multiplication for all other terms and sum them up if there are non-zero coefficients
        auto qcopy = divqrVec[0]->q;
        qcopy.resize(k);

        if (Degree(qcopy) > 0) {
            std::vector<Ciphertext<DCRTPoly>> ctxs(Degree(qcopy));
            std::vector<std::vector<double>> weights(divqrVec.size());

            for (size_t j = 0; j < divqrVec.size(); ++j) {
                std::vector<double> weightsForLevel;
                for (uint32_t i = 0; i < Degree(qcopy); ++i) {
                    ctxs[i] = T[i];
                    weightsForLevel.push_back(divqrVec[j]->q[i + 1]);
                }
                weights[j] = weightsForLevel;
            }

            qu = internalEvalLinearWSumMutableBatch(ctxs, weights);

            // the highest order coefficient will always be a power of two up to 2^{m-1} because q is "monic" but the Chebyshev rule adds a factor of 2
            // we don't need to increase the depth by multiplying the highest order coefficient, but instead checking and summing, since we work with m <= 4.
            Ciphertext<DCRTPoly> sum = T[k - 1]->Clone();
            uint32_t limit = log2(ToReal(divqrVec[0]->q.back()));
            for (uint32_t i = 0; i < limit; ++i) {
                sum = cc->EvalAdd(sum, sum);
            }

            qu = cc->EvalAdd(qu, sum);
        }
        else {
            Ciphertext<DCRTPoly> sum = T[k - 1]->Clone();
            uint32_t limit = log2(ToReal(divqrVec[0]->q.back()));
            for (uint32_t i = 0; i < limit; ++i) {
                sum = cc->EvalAdd(sum, sum);
            }

            qu = sum;
        }

        // adds the free term (at x^0)
        std::vector<double> freeTerms;
        for (size_t i = 0; i < divqrVec.size(); i++) {
            freeTerms.push_back(divqrVec[i]->q.front() / 2.0);
        }

        Plaintext ptxt = cc->MakeCKKSPackedPlaintext(freeTerms, 1, qu->GetLevel(),
                                                     nullptr, qu->GetSlots());
        qu = cc->EvalAdd(qu, ptxt);
        // The number of levels of qu is the same as the number of levels of T[k-1] + 1.
        // Will only get here when m = 2, so the number of levels of qu and T2[m-1] will be the same.
    }

    Ciphertext<DCRTPoly> su;

    if (Degree(s2Vec[0]) > k) {
        su = InnerEvalChebyshevPSBatch(T[0], s2Vec, k, m - 1, T, T2);
    }
    else {
        // ds = k from construction
        // perform scalar multiplication for all other terms and sum them up if there are non-zero coefficients
        auto scopy = s2Vec[0];
        scopy.resize(k);
        if (Degree(scopy) > 0) {
            std::vector<Ciphertext<DCRTPoly>> ctxs(Degree(scopy));
            std::vector<std::vector<double>> weights(s2Vec.size());

            for (uint32_t j = 0; j < s2Vec.size(); j++) {
                std::vector<double> weightsForLevel;
                for (uint32_t i = 0; i < Degree(scopy); ++i) {
                    ctxs[i] = T[i];
                    weightsForLevel.push_back(s2Vec[j][i + 1]);
                }
                weights[j] = weightsForLevel;
            }

            su = cc->EvalLinearWSumMutableBatch(ctxs, weights);

            // the highest order coefficient will always be 1 because s2 is monic.
            cc->EvalAddInPlace(su, T[k - 1]);
        }
        else {
            su = T[k - 1];
        }

        // adds the free term (at x^0)
        std::vector<double> freeTerms;
        for (uint32_t i = 0; i < s2Vec.size(); i++) {
            freeTerms.push_back(s2Vec[i].front() / 2.0);
        }

        Plaintext ptxt = cc->MakeCKKSPackedPlaintext(freeTerms, 1, su->GetLevel(), nullptr, su->GetSlots());
        cc->EvalAddInPlace(su, ptxt);

        // The number of levels of su is the same as the number of levels of T[k-1] + 1.
        // Will only get here when m = 2, so need to reduce the number of levels by 1.
    }

    // TODO : Andrey : here is different from 895 line
    // Reduce number of levels of su to number of levels of T2km1.
    //  cc->LevelReduceInPlace(su, nullptr);

    Ciphertext<DCRTPoly> result;

    if (flag_c) {
        result = cc->EvalAdd(T2[m - 1], cu);
    }
    else {
        std::vector<double> terms;
        for (size_t i = 0; i < divcsVec.size(); i++) {
            terms.push_back(divcsVec[i]->q.front() / 2.0);
        }

        Plaintext ptxt = cc->MakeCKKSPackedPlaintext(terms, 1, T2[m - 1]->GetLevel(), nullptr, T2[m - 1]->GetSlots());
        result = cc->EvalAdd(T2[m - 1], ptxt);
    }

    result = cc->EvalMult(result, qu);
    cc->ModReduceInPlace(result);

    cc->EvalAddInPlace(result, su);
    cc->EvalSubInPlace(result, T2km1);

    return result;
}

template <typename VectorDataType>
static inline Ciphertext<DCRTPoly> internalEvalChebyshevSeriesPSBatchRepeatedWithPrecomp(
    std::shared_ptr<seriesPowers<DCRTPoly>> ctxtPolys,
    const std::vector<std::vector<VectorDataType>>& batchOfCoefficients, int repetitions) {

    std::vector<std::vector<double>> batchOfCoefficientsDouble;
    batchOfCoefficientsDouble.reserve(batchOfCoefficients.size());

    // This is required as throughout the function we will create CKKS plaintexts,
    // but the encoding function requires doubles (and not complex values, accepted by
    // VectorDataType)
    if constexpr (std::is_same_v<VectorDataType, double>) {
        for (const auto& innerVec : batchOfCoefficients) {
            std::vector<double> converted;
            converted.reserve(innerVec.size());

            for (const auto& value : innerVec) {
                converted.push_back(static_cast<double>(value));
            }
            batchOfCoefficientsDouble.push_back(std::move(converted));
        }
    } else {
        OPENFHE_THROW(std::string("Not implemented for type: ") + typeid(VectorDataType).name());
    }


    std::vector<std::shared_ptr<longDiv<double>>> divcsVec;
    std::vector<std::shared_ptr<longDiv<double>>> divqrVec;
    std::vector<std::vector<double>> s2Vec;

    /*
     * Check on the batch of coefficients: they must all be of the same size
     */
    std::size_t size = batchOfCoefficients[0].size();
    for (const auto& inner : batchOfCoefficients) {
        if (inner.size() != size) {
            OPENFHE_THROW("The batch must contain coefficients of the same size!");
        }
    }

    auto T     = ctxtPolys->powersRe;
    auto T2    = ctxtPolys->powers2Re;
    auto T2km1 = ctxtPolys->power2km1Re;
    auto k     = ctxtPolys->k;
    auto m     = ctxtPolys->m;

    repetitions = T.front()->GetSlots() / batchOfCoefficientsDouble.size();

    // Compute k*2^{m-1}-k because we use it a lot
    uint32_t k2m2k = k * (1 << (m - 1)) - k;

    for (size_t b = 0; b < batchOfCoefficientsDouble.size(); b++) {
        auto f2 = batchOfCoefficientsDouble[b];
        auto n  = Degree(f2);
        f2.resize(n + 1);

        // Add T^{k(2^m - 1)}(y) to the polynomial that has to be evaluated
        f2.resize(2 * k2m2k + k + 1, 0.0);
        f2.back() = 1;

        // Divide f2 by T^{k*2^{m-1}}
        std::vector<double> Tkm(k2m2k + k + 1);
        Tkm.back() = 1;

        auto divqr = LongDivisionChebyshev(f2, Tkm);

        // Subtract x^{k(2^{m-1} - 1)} from r
        auto r2 = divqr->r;
        if (static_cast<int32_t>(k2m2k - Degree(r2)) <= 0) {
            r2[static_cast<int32_t>(k2m2k)] -= 1;
            r2.resize(Degree(r2) + 1);
        }
        else {
            r2.resize(static_cast<int32_t>(k2m2k + 1));
            r2.back() = -1;
        }

        // Divide r2 by q
        auto divcs = LongDivisionChebyshev(r2, divqr->q);

        // Add x^{k(2^{m-1} - 1)} to s
        auto s2 = divcs->r;
        s2.resize(k2m2k + 1);
        s2.back() = 1;

        divqrVec.push_back(divqr);
        divcsVec.push_back(divcs);
        s2Vec.push_back(s2);
    }

    auto cc = T[0]->GetCryptoContext();

    // Evaluate c at u
    Ciphertext<DCRTPoly> cu;

    // We use the degree of the first one, they should all be the same anyways
    uint32_t dc = Degree(divcsVec[0]->q);
    bool flag_c = false;
    if (dc >= 1) {
        if (dc == 1) {
            // if (IsNotEqualOne(divcs->q[1])) {
            // NOTE: We can not optimize anymore since it is not a constant, so we always perform the product
            std::vector<double> coeffs;
            for (size_t i = 0; i < divcsVec.size(); i++) coeffs.push_back(divcsVec[i]->q[1]);

            //std::cout << coeffs.size() << std::endl;

            std::vector<double> coeffsRep;
            for (int i = 0; i < repetitions; ++i) {
                coeffsRep.insert(coeffsRep.end(), coeffs.begin(), coeffs.end());
            }

            cu = cc->EvalMult(T.front(), cc->MakeCKKSPackedPlaintext(coeffsRep, 1, T.front()->GetLevel(),
                                                                     nullptr, T.front()->GetSlots()));
            cc->ModReduceInPlace(cu);
        }
        else {
            std::vector<std::vector<double>> batchOfWeights;
            std::vector<Ciphertext<DCRTPoly>> ctxs(dc);

            for (size_t i = 0; i < dc; i++) {
                ctxs[i] = T[i];
            }

            for (size_t j = 0; j < divcsVec.size(); j++) {
                std::vector<double> weights(dc);
                for (size_t i = 0; i < dc; i++) {
                    weights[i] = divcsVec[j]->q[i + 1];
                }

                batchOfWeights.push_back(weights);
            }

            cu = internalEvalLinearWSumMutableBatch(ctxs, batchOfWeights);
        }

        // adds the free term (at x^0)
        std::vector<double> freeTerms;
        for (size_t i = 0; i < divcsVec.size(); i++) {
            freeTerms.push_back(divcsVec[i]->q.front() / 2.0);
        }

        std::vector<double> freeTermsRep;
        for (int i = 0; i < repetitions; ++i) {
            freeTermsRep.insert(freeTermsRep.end(), freeTerms.begin(), freeTerms.end());
        }

        Plaintext freeTermsPtxt = cc->MakeCKKSPackedPlaintext(freeTermsRep, 1, cu->GetLevel(),
                                                              nullptr, cu->GetSlots());
        cu = cc->EvalAdd(cu, freeTermsPtxt);

        // TODO : Andrey why not T2[m-1]->GetLevel() instead?
        // Need to reduce levels to the level of T2[m-1].
        //    uint32_t levelDiff = y->GetLevel() - cu->GetLevel() + ceil(log2(k)) + m - 1;
        //    cc->LevelReduceInPlace(cu, nullptr, levelDiff);

        flag_c = true;
    }

    // Evaluate q and s2 at u. If their degrees are larger than k, then recursively apply the Paterson-Stockmeyer algorithm.
    Ciphertext<DCRTPoly> qu;

    // Again, degrees should all be the same for the different divqrVec
    if (Degree(divqrVec[0]->q) > k) {
        std::vector<std::vector<double>> coeffs;
        for (size_t i = 0; i < divqrVec.size(); i++) {
            coeffs.push_back(divqrVec[i]->q);
        }
        qu = InnerEvalChebyshevPSBatchRepeated(T[0], coeffs, k, m - 1, T, T2);
    }
    else {
        // dq = k from construction
        // perform scalar multiplication for all other terms and sum them up if there are non-zero coefficients
        auto qcopy = divqrVec[0]->q;
        qcopy.resize(k);

        if (Degree(qcopy) > 0) {
            std::vector<Ciphertext<DCRTPoly>> ctxs(Degree(qcopy));
            std::vector<std::vector<double>> weights(divqrVec.size());

            for (size_t j = 0; j < divqrVec.size(); ++j) {
                std::vector<double> weightsForLevel;
                for (uint32_t i = 0; i < Degree(qcopy); ++i) {
                    ctxs[i] = T[i];
                    weightsForLevel.push_back(divqrVec[j]->q[i + 1]);
                }
                weights[j] = weightsForLevel;
            }

            qu = internalEvalLinearWSumMutableBatch(ctxs, weights);

            // the highest order coefficient will always be a power of two up to 2^{m-1} because q is "monic" but the Chebyshev rule adds a factor of 2
            // we don't need to increase the depth by multiplying the highest order coefficient, but instead checking and summing, since we work with m <= 4.
            Ciphertext<DCRTPoly> sum = T[k - 1]->Clone();
            uint32_t limit = log2(ToReal(divqrVec[0]->q.back()));
            for (uint32_t i = 0; i < limit; ++i) {
                sum = cc->EvalAdd(sum, sum);
            }

            qu = cc->EvalAdd(qu, sum);
        }
        else {
            Ciphertext<DCRTPoly> sum = T[k - 1]->Clone();
            uint32_t limit = log2(ToReal(divqrVec[0]->q.back()));
            for (uint32_t i = 0; i < limit; ++i) {
                sum = cc->EvalAdd(sum, sum);
            }

            qu = sum;
        }

        // adds the free term (at x^0)
        std::vector<double> freeTerms;
        for (size_t i = 0; i < divqrVec.size(); i++) {
            freeTerms.push_back(divqrVec[i]->q.front() / 2.0);
        }

        std::vector<double> freeTermsRep;
        for (int i = 0; i < repetitions; ++i) {
            freeTermsRep.insert(freeTermsRep.end(), freeTerms.begin(), freeTerms.end());
        }

        Plaintext ptxt = cc->MakeCKKSPackedPlaintext(freeTermsRep, 1, qu->GetLevel(),
                                                     nullptr, qu->GetSlots());
        qu = cc->EvalAdd(qu, ptxt);
        // The number of levels of qu is the same as the number of levels of T[k-1] + 1.
        // Will only get here when m = 2, so the number of levels of qu and T2[m-1] will be the same.
    }

    Ciphertext<DCRTPoly> su;

    if (Degree(s2Vec[0]) > k) {
        su = InnerEvalChebyshevPSBatchRepeated(T[0], s2Vec, k, m - 1, T, T2);
    }
    else {
        // ds = k from construction
        // perform scalar multiplication for all other terms and sum them up if there are non-zero coefficients
        auto scopy = s2Vec[0];
        scopy.resize(k);
        if (Degree(scopy) > 0) {
            std::vector<Ciphertext<DCRTPoly>> ctxs(Degree(scopy));
            std::vector<std::vector<double>> weights(s2Vec.size());

            for (uint32_t j = 0; j < s2Vec.size(); j++) {
                std::vector<double> weightsForLevel;
                for (uint32_t i = 0; i < Degree(scopy); ++i) {
                    ctxs[i] = T[i];
                    weightsForLevel.push_back(s2Vec[j][i + 1]);
                }
                weights[j] = weightsForLevel;
            }

            su = cc->EvalLinearWSumMutableBatch(ctxs, weights);

            // the highest order coefficient will always be 1 because s2 is monic.
            cc->EvalAddInPlace(su, T[k - 1]);
        }
        else {
            su = T[k - 1];
        }

        // adds the free term (at x^0)
        std::vector<double> freeTerms;
        for (uint32_t i = 0; i < s2Vec.size(); i++) {
            freeTerms.push_back(s2Vec[i].front() / 2.0);
        }

        std::vector<double> freeTermsRep;
        for (int i = 0; i < repetitions; ++i) {
            freeTermsRep.insert(freeTermsRep.end(), freeTerms.begin(), freeTerms.end());
        }

        Plaintext ptxt = cc->MakeCKKSPackedPlaintext(freeTermsRep, 1, su->GetLevel(), nullptr, su->GetSlots());
        cc->EvalAddInPlace(su, ptxt);

        // The number of levels of su is the same as the number of levels of T[k-1] + 1.
        // Will only get here when m = 2, so need to reduce the number of levels by 1.
    }

    // TODO : Andrey : here is different from 895 line
    // Reduce number of levels of su to number of levels of T2km1.
    //  cc->LevelReduceInPlace(su, nullptr);

    Ciphertext<DCRTPoly> result;

    if (flag_c) {
        result = cc->EvalAdd(T2[m - 1], cu);
    }
    else {
        std::vector<double> terms;
        for (size_t i = 0; i < divcsVec.size(); i++) {
            terms.push_back(divcsVec[i]->q.front() / 2.0);
        }

        std::vector<double> termsRep;
        for (int i = 0; i < repetitions; ++i) {
            termsRep.insert(termsRep.end(), terms.begin(), terms.end());
        }

        Plaintext ptxt = cc->MakeCKKSPackedPlaintext(termsRep, 1, T2[m - 1]->GetLevel(), nullptr, T2[m - 1]->GetSlots());
        result = cc->EvalAdd(T2[m - 1], ptxt);
    }

    result = cc->EvalMult(result, qu);
    cc->ModReduceInPlace(result);

    cc->EvalAddInPlace(result, su);
    cc->EvalSubInPlace(result, T2km1);

    return result;
}

std::shared_ptr<seriesPowers<DCRTPoly>> AdvancedSHECKKSRNS::EvalChebyPolys(ConstCiphertext<DCRTPoly>& x,
                                                                           const std::vector<int64_t>& coefficients,
                                                                           double a, double b) const {
    uint32_t d = Degree(coefficients);
    return (d < 5) ? internalEvalChebyPolysLinear(x, coefficients, a, b) : internalEvalChebyPolysPS(x, d, a, b);
}
std::shared_ptr<seriesPowers<DCRTPoly>> AdvancedSHECKKSRNS::EvalChebyPolys(ConstCiphertext<DCRTPoly>& x,
                                                                           const std::vector<double>& coefficients,
                                                                           double a, double b) const {
    uint32_t d = Degree(coefficients);
    return (d < 5) ? internalEvalChebyPolysLinear(x, coefficients, a, b) : internalEvalChebyPolysPS(x, d, a, b);
}
std::shared_ptr<seriesPowers<DCRTPoly>> AdvancedSHECKKSRNS::EvalChebyPolys(
    ConstCiphertext<DCRTPoly>& x, const std::vector<std::complex<double>>& coefficients, double a, double b) const {
    uint32_t d = Degree(coefficients);
    return (d < 5) ? internalEvalChebyPolysLinear(x, coefficients, a, b) : internalEvalChebyPolysPS(x, d, a, b);
}

Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalChebyshevSeries(ConstCiphertext<DCRTPoly>& x,
                                                             const std::vector<int64_t>& coeffs, double a,
                                                             double b) const {
    return (Degree(coeffs) < 5) ? EvalChebyshevSeriesLinear(x, coeffs, a, b) : EvalChebyshevSeriesPS(x, coeffs, a, b);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalChebyshevSeries(ConstCiphertext<DCRTPoly>& x,
                                                             const std::vector<double>& coeffs, double a,
                                                             double b) const {
    return (Degree(coeffs) < 5) ? EvalChebyshevSeriesLinear(x, coeffs, a, b) : EvalChebyshevSeriesPS(x, coeffs, a, b);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalChebyshevSeries(ConstCiphertext<DCRTPoly>& x,
                                                             const std::vector<std::complex<double>>& coeffs, double a,
                                                             double b) const {
    return (Degree(coeffs) < 5) ? EvalChebyshevSeriesLinear(x, coeffs, a, b) : EvalChebyshevSeriesPS(x, coeffs, a, b);
}

Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalChebyshevSeriesWithPrecomp(
    std::shared_ptr<seriesPowers<DCRTPoly>> ctxtPowers, const std::vector<int64_t>& coeffs) const {
    return (Degree(coeffs) < 5) ? internalEvalChebyshevSeriesLinearWithPrecomp(ctxtPowers->powersRe, coeffs) :
                                  internalEvalChebyshevSeriesPSWithPrecomp(ctxtPowers, coeffs);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalChebyshevSeriesWithPrecomp(
    std::shared_ptr<seriesPowers<DCRTPoly>> ctxtPowers, const std::vector<double>& coeffs) const {
    return (Degree(coeffs) < 5) ? internalEvalChebyshevSeriesLinearWithPrecomp(ctxtPowers->powersRe, coeffs) :
                                  internalEvalChebyshevSeriesPSWithPrecomp(ctxtPowers, coeffs);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalChebyshevSeriesWithPrecomp(
    std::shared_ptr<seriesPowers<DCRTPoly>> ctxtPowers, const std::vector<std::complex<double>>& coeffs) const {
    return (Degree(coeffs) < 5) ? internalEvalChebyshevSeriesLinearWithPrecomp(ctxtPowers->powersRe, coeffs) :
                                  internalEvalChebyshevSeriesPSWithPrecomp(ctxtPowers, coeffs);
}

Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalChebyshevSeriesLinear(ConstCiphertext<DCRTPoly>& x,
                                                                   const std::vector<int64_t>& coeffs, double a,
                                                                   double b) const {
    return internalEvalChebyshevSeriesLinearWithPrecomp(internalEvalChebyPolysLinear(x, coeffs, a, b)->powersRe,
                                                        coeffs);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalChebyshevSeriesLinear(ConstCiphertext<DCRTPoly>& x,
                                                                   const std::vector<double>& coeffs, double a,
                                                                   double b) const {
    return internalEvalChebyshevSeriesLinearWithPrecomp(internalEvalChebyPolysLinear(x, coeffs, a, b)->powersRe,
                                                        coeffs);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalChebyshevSeriesLinear(ConstCiphertext<DCRTPoly>& x,
                                                                   const std::vector<std::complex<double>>& coeffs,
                                                                   double a, double b) const {
    return internalEvalChebyshevSeriesLinearWithPrecomp(internalEvalChebyPolysLinear(x, coeffs, a, b)->powersRe,
                                                        coeffs);
}

Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalChebyshevSeriesPS(ConstCiphertext<DCRTPoly>& x,
                                                               const std::vector<int64_t>& coeffs, double a,
                                                               double b) const {
    return internalEvalChebyshevSeriesPSWithPrecomp(internalEvalChebyPolysPS(x, Degree(coeffs), a, b), coeffs);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalChebyshevSeriesPS(ConstCiphertext<DCRTPoly>& x,
                                                               const std::vector<double>& coeffs, double a,
                                                               double b) const {
    return internalEvalChebyshevSeriesPSWithPrecomp(internalEvalChebyPolysPS(x, Degree(coeffs), a, b), coeffs);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalChebyshevSeriesPS(ConstCiphertext<DCRTPoly>& x,
                                                               const std::vector<std::complex<double>>& coeffs,
                                                               double a, double b) const {
    return internalEvalChebyshevSeriesPSWithPrecomp(internalEvalChebyPolysPS(x, Degree(coeffs), a, b), coeffs);
}

//For batched SIMD Chebyshev
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalChebyshevSeriesPSBatch(ConstCiphertext<DCRTPoly>& x,
                                                               const std::vector<std::vector<int64_t>>& batchOfCoefficients,
                                                               double a, double b) const {
    if (batchOfCoefficients.size() != x->GetSlots())
        OPENFHE_THROW("The set of coefficients must be as large as the number of slots of the input ciphertext");

    auto deg = Degree(batchOfCoefficients[0]);
    for (int i = 0; i < batchOfCoefficients.size(); i++) {
        if (Degree(batchOfCoefficients[i]) != deg) {
            OPENFHE_THROW("The polynomials must have all the same degrees");
        }
    }

    return internalEvalChebyshevSeriesPSBatchWithPrecomp(internalEvalChebyPolysPS(x, Degree(batchOfCoefficients[0]), a, b), batchOfCoefficients);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalChebyshevSeriesPSBatch(ConstCiphertext<DCRTPoly>& x,
                                                                      const std::vector<std::vector<double>>& batchOfCoefficients,
                                                                      double a, double b) const {
    if (batchOfCoefficients.size() != x->GetSlots())
        OPENFHE_THROW("The set of coefficients must be as large as the number of slots of the input ciphertext");

    auto deg = Degree(batchOfCoefficients[0]);
    for (int i = 0; i < batchOfCoefficients.size(); i++) {
        if (Degree(batchOfCoefficients[i]) != deg) {
            OPENFHE_THROW("The polynomials must have all the same degrees");
        }
    }

    return internalEvalChebyshevSeriesPSBatchWithPrecomp(internalEvalChebyPolysPS(x, Degree(batchOfCoefficients[0]), a, b), batchOfCoefficients);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalChebyshevSeriesPSBatch(
    ConstCiphertext<DCRTPoly>& x, const std::vector<std::vector<std::complex<double>>>& batchOfCoefficients, double a, double b)
    const {
    if (batchOfCoefficients.size() != x->GetSlots())
        OPENFHE_THROW("The set of coefficients must be as large as the number of slots of the input ciphertext");

    auto deg = Degree(batchOfCoefficients[0]);
    for (int i = 0; i < batchOfCoefficients.size(); i++) {
        if (Degree(batchOfCoefficients[i]) != deg) {
            OPENFHE_THROW("The polynomials must have all the same degrees");
        }
    }

    return internalEvalChebyshevSeriesPSBatchWithPrecomp(internalEvalChebyPolysPS(x, Degree(batchOfCoefficients[0]), a, b), batchOfCoefficients);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalChebyshevSeriesPSBatchRepeated(ConstCiphertext<DCRTPoly>& x,
                                                                    const std::vector<std::vector<int64_t>>& batchOfCoefficients,
                                                                    double a, double b, int repetitions) const {
    return internalEvalChebyshevSeriesPSBatchRepeatedWithPrecomp(internalEvalChebyPolysPS(x, Degree(batchOfCoefficients[0]), a, b), batchOfCoefficients, repetitions);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalChebyshevSeriesPSBatchRepeated(ConstCiphertext<DCRTPoly>& x,
                                                                            const std::vector<std::vector<std::complex<double>>>& batchOfCoefficients,
                                                                            double a, double b, int repetitions) const {
    return internalEvalChebyshevSeriesPSBatchRepeatedWithPrecomp(internalEvalChebyPolysPS(x, Degree(batchOfCoefficients[0]), a, b), batchOfCoefficients, repetitions);
}
Ciphertext<DCRTPoly> AdvancedSHECKKSRNS::EvalChebyshevSeriesPSBatchRepeated(ConstCiphertext<DCRTPoly>& x,
                                                                            const std::vector<std::vector<double>>& batchOfCoefficients,
                                                                            double a, double b, int repetitions) const {
    return internalEvalChebyshevSeriesPSBatchRepeatedWithPrecomp(internalEvalChebyPolysPS(x, Degree(batchOfCoefficients[0]), a, b), batchOfCoefficients, repetitions);
}


//------------------------------------------------------------------------------
// EVAL LINEAR TRANSFORMATION
//------------------------------------------------------------------------------

}  // namespace lbcrypto
