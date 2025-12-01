#include <iostream>
#include "openfhe.h"

using namespace lbcrypto;
using namespace std;

// Define aliases for clarity
using Ptxt = Plaintext;
using Ctxt = Ciphertext<DCRTPoly>;

int main(int argc, char *argv[]) {
    int num_slots = 2;

    CCParams<CryptoContextCKKSRNS> parameters;

    parameters.SetSecretKeyDist(lbcrypto::UNIFORM_TERNARY);

    int dcrtBits = 40;
    int firstMod = 45;

    parameters.SetSecurityLevel(lbcrypto::HEStd_NotSet);
    parameters.SetRingDim(1 << 12);

    cout << "N: " << parameters.GetRingDim() << endl << endl;

    parameters.SetBatchSize(num_slots);

    ScalingTechnique rescaleTech = FLEXIBLEAUTO;

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetFirstModSize(firstMod);
    parameters.SetMultiplicativeDepth(12);


    CryptoContext<DCRTPoly> context = GenCryptoContext(parameters);
    context->Enable(PKE);
    context->Enable(KEYSWITCH);
    context->Enable(LEVELEDSHE);
    context->Enable(ADVANCEDSHE);

    KeyPair<DCRTPoly> key_pair = context->KeyGen();
    context->EvalMultKeyGen(key_pair.secretKey);

    vector<double> input = {0.1, 0.2};

    Ptxt p = context->MakeCKKSPackedPlaintext(input, 1, 0, nullptr, num_slots);
    Ctxt c = context->Encrypt(p, key_pair.publicKey);

    vector<double> coeffs1 = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

    cout << "Test with 6 coefficients " << coeffs1 << " over x = 0.1, 0.2" << endl;

    Ctxt cres = context->EvalChebyshevSeriesPS(c, coeffs1, -1, 1);

    Ptxt res;
    context->Decrypt(cres, key_pair.secretKey, &res);
    cout << res << endl;

    vector<double> coeffs2 = {4.0, 1.5, 2.0, 0.5, 0.1, 2.0};

    cout << "Test with 6 coefficients " << coeffs2 << " over x = 0.1, 0.2" << endl;

    cres = context->EvalChebyshevSeriesPS(c, coeffs2, -1, 1);
    context->Decrypt(cres, key_pair.secretKey, &res);
    cout << res << endl;

    cout << "Now parallelizing the two polynomials in two slots" << endl;

    cres = context->EvalChebyshevSeriesPSBatch(c, {coeffs1, coeffs2}, -1, 1);
    context->Decrypt(cres, key_pair.secretKey, &res);
    cout << res << endl;

    /*
     * TEST WITH MANY COEFFS
     */
    coeffs1 = {3.4, 1.2, 0.8, 1.3, 0.0, 2.3, 1.1, 0.6, 0.3, 0.6, 3.0, 2.1, 0.3, 2.7, 2.9, 1.7, 2.2, 0.5, 0.5, 0.6, 1.6, 0.8, 0.6, 0.8, 1.8, 1.4, 1.6, 1.7, 2.5, 2.2, 0.2, 1.6, 1.9, 1.2, 1.9, 0.1, 1.0, 0.1, 1.4, 2.5, 2.9, 2.9, 2.7, 0.2, 0.0, 0.2, 2.0, 0.8, 1.3, 2.0};
    cout << "Test with 50 coefficients " << coeffs1 << " over x = 0.1, 0.2" << endl;

    cres = context->EvalChebyshevSeriesPS(c, coeffs1, -1, 1);

    context->Decrypt(cres, key_pair.secretKey, &res);
    cout << res << endl;

    coeffs2 = {-0.5, 1.8, 1.3, 1.9, -1.1, 1.7, -2.0, 1.9, -1.8, -1.8, 1.8, -0.4, -0.5, -1.1, 1.5, -1.0, -0.3, 0.4, -1.2, -0.4, -0.2, -0.2, -1.8, 1.6, -1.4, -1.4, 0.1, -1.3, -1.7, -0.7, -0.6, 0.7, 1.8, -0.5, -0.3, -1.8, -2.0, -0.7, 1.0, 0.2, 1.5, 0.0, -1.1, 1.2, -0.4, 2.0, 0.5, -2.0, -1.6, 0.7};
    cout << "Test with 50 coefficients " << coeffs2 << " over x = 0.1, 0.2" << endl;

    cres = context->EvalChebyshevSeriesPS(c, coeffs2, -1, 1);
    context->Decrypt(cres, key_pair.secretKey, &res);
    cout << res << endl;

    cout << "Now parallelizing the two polynomials in two slots" << endl;

    cres = context->EvalChebyshevSeriesPSBatch(c, {coeffs1, coeffs2}, -1, 1);
    context->Decrypt(cres, key_pair.secretKey, &res);
    cout << res << endl;


}

