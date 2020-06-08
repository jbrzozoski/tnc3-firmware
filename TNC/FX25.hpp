#pragma once

// #include "HdlcFrame.hpp"
// #include "main.h"
// #include <cmsis_os.h>
// #include <cstdint>
#include <cstring>
#include <cstdio>

//namespace mobilinkd { namespace tnc {

// The ReedSolomon engine takes several constant configuration values.
// Originally I was going to make them parameters to the constructor, but many of them
// affect the sizes of the internal data structures needed.
// It is wildly improbably that someone might want to try a different parameter value at runtime, so
// instead of playing games with allocating memory as part of the constructor, I just moved all the
// parameters into the template.  I added static asserts for validity checking of the parameters.
template <typename DTYPE, uint16_t SYMSIZE, uint16_t GFPOLY, uint8_t FCR, uint8_t PRIM, uint16_t NROOTS>
struct rs_config {
    // DTYPE: data type being operated on (expected to be uint8_t)
    // SYMSIZE: Bits per symbol (aka MM in some reference code) (expected to be 8)
    // GFPOLY: Field generator polynomial coefficients
    // FCR: First consecutive root, index form
    // PRIM: Primitive element, index form
    // NROOTS: Number of generator roots = number of parity symbols
    static_assert(SYMSIZE <= (8 * sizeof(DTYPE)));
    static_assert(FCR < (1 << SYMSIZE));
    static_assert((PRIM != 0) && (PRIM < (1 << SYMSIZE)));
    static_assert(NROOTS < (1 << SYMSIZE));
    const uint16_t NN = ((1<<SYMSIZE)-1); /* Symbols per block */
    const uint16_t A0 = NN;
    uint8_t alpha_to[(sizeof(DTYPE) * (((1<<SYMSIZE)-1) + 1))]; /* log lookup table */
    uint8_t index_of[(sizeof(DTYPE) * (((1<<SYMSIZE)-1) + 1))]; /* Antilog lookup table */
    uint8_t genpoly[(sizeof(DTYPE) * (NROOTS + 1))];    /* Generator polynomial */
    uint8_t iprim;       /* prim-th root of 1, index form */

    inline uint16_t modnn(uint16_t x) const {
        while (x >= NN) {
            x -= NN;
            x = (x >> SYMSIZE) + (x & NN);
        }
        return x;
    }

    rs_config()
    : alpha_to(), index_of(), genpoly(), iprim()
    {
        int i, j, sr, root;

        /* Generate Galois field lookup tables */
        index_of[0] = A0; /* log(zero) = -inf */
        alpha_to[A0] = 0; /* alpha**-inf = 0 */
        sr = 1;
        for (i = 0; i < NN; i++) {
            index_of[sr] = i;
            alpha_to[i] = sr;
            sr <<= 1;
            if (sr & (1 << SYMSIZE)) {
                sr ^= GFPOLY;
            }
            sr &= NN;
        }
        if (sr != 1) {
            /* field generator polynomial is not primitive! */
            // Should find a better way ...???
            // I don't think a better exception class would help (std::domain_error is likely if we do, though...)
            throw 0;
        }

        /* Form RS code generator polynomial from its roots. */
        /* Find prim-th root of 1, used in decoding */
        iprim = 1;
        while ((iprim % PRIM) != 0) {
            iprim += NN;
        }
        iprim = iprim / PRIM;
        genpoly[0] = 1;
        for (i = 0, root = FCR * PRIM; i < NROOTS; i++, root += PRIM) {
            genpoly[i + 1] = 1;

            /* Multiply genpoly[] by  @**(root + x) */
            for (j = i; j > 0; j--) {
                if (genpoly[j] != 0) {
                    genpoly[j] = genpoly[j - 1] ^ alpha_to[modnn(index_of[genpoly[j]] + root)];
                } else {
                    genpoly[j] = genpoly[j - 1];
                }
            }
            /* genpoly[0] can never be zero */
            genpoly[0] = alpha_to[modnn(index_of[genpoly[0]] + root)];
        }

        /* convert genpoly[] to index form for quicker encoding */
        for (i = 0; i <= NROOTS; i++) {
            genpoly[i] = index_of[genpoly[i]];
        }

        // diagnostic prints
        printf("Alpha To:\n");
        for (i = 0; i < sizeof(DTYPE) * (NN + 1); i++) {
            printf("%02x ", alpha_to[i]);
            if ((i % 16) == 15) {
                printf("\n");
            }
        }
        printf("\n");

        printf("Index Of:\n");
        for (i = 0; i < sizeof(DTYPE) * (NN + 1); i++) {
            printf("%02x ", index_of[i]);
            if ((i % 16) == 15) {
                printf("\n");
            }
        }
        printf("\n");

        printf("GenPoly:\n");
        for (i = 0; i <= NROOTS; i++) {
            printf("%02x ", genpoly[i]);
            if ((i % 16) == 15) {
                printf("\n");
            }
        }
        printf("\n");
    }

    void ENCODE_RS(DTYPE * data, DTYPE * bb)
    {
        int i, j;
        DTYPE feedback;

        // clear out the FEC data area
        memset(bb, 0, NROOTS * sizeof(DTYPE));

        for(i = 0; i < (NN-NROOTS); i++) {
            feedback = index_of[data[i] ^ bb[0]];
            /* feedback term is non-zero */
            if(feedback != A0) {
                for(j = 1; j < NROOTS; j++) {
                    bb[j] ^= alpha_to[modnn(feedback + genpoly[NROOTS-j])];
                }
            }
            /* Shift */
            memmove(&bb[0], &bb[1], sizeof(DTYPE)*(NROOTS-1));
            if(feedback != A0) {
                bb[NROOTS-1] = alpha_to[modnn(feedback + genpoly[0])];
            } else {
                bb[NROOTS-1] = 0;
            }
        }
    }
};

// struct fx25 {
//     struct rs_config<uint8_t, 8, 0x11d, 1, 1, 16> RS_255_239;
//     struct rs_config<uint8_t, 8, 0x11d, 1, 1, 32> RS_255_223;
//     struct rs_config<uint8_t, 8, 0x11d, 1, 1, 64> RS_255_191;
//
//     fx25() {
//         // Verify integrity of tables and assumptions.
//         // This also does a quick check for the popcount function.
//         for (int j = 0; j < 16; j++) {
//             for (int k = 0; k < 16; k++) {
//                 if (j == k) {
//                     assert(__builtin_popcountll(tags[j].value ^ tags[k].value) == 0);
//                 } else {
//                     assert(__builtin_popcountll(tags[j].value ^ tags[k].value) == 32);
//                 }
//             }
//         }
//
//         for (int j = CTAG_MIN; j <= CTAG_MAX; j++) {
//             assert(tags[j].n_block_radio - tags[j].k_data_radio == Tab[tags[j].itab].NROOTS);
//             assert(tags[j].n_block_rs - tags[j].k_data_rs == Tab[tags[j].itab].NROOTS);
//             assert(tags[j].n_block_rs == FX25_BLOCK_SIZE);
//         }
//
//         assert(fx25_pick_mode(1, 239) == 1);
//         assert(fx25_pick_mode(1, 240) == -1);
//
//         assert(fx25_pick_mode(5, 223) == 5);
//         assert(fx25_pick_mode(5, 224) == -1);
//
//         assert(fx25_pick_mode(9, 191) == 9);
//         assert(fx25_pick_mode(9, 192) == -1);
//
//         assert(fx25_pick_mode(16, 32) == 4);
//         assert(fx25_pick_mode(16, 64) == 3);
//         assert(fx25_pick_mode(16, 128) == 2);
//         assert(fx25_pick_mode(16, 239) == 1);
//         assert(fx25_pick_mode(16, 240) == -1);
//
//         assert(fx25_pick_mode(32, 32) == 8);
//         assert(fx25_pick_mode(32, 64) == 7);
//         assert(fx25_pick_mode(32, 128) == 6);
//         assert(fx25_pick_mode(32, 223) == 5);
//         assert(fx25_pick_mode(32, 234) == -1);
//
//         assert(fx25_pick_mode(64, 64) == 11);
//         assert(fx25_pick_mode(64, 128) == 10);
//         assert(fx25_pick_mode(64, 191) == 9);
//         assert(fx25_pick_mode(64, 192) == -1);
//     }
//
//     void encode(IoFrame *frame) {
//     }
// }
//
// struct correlation_tag_s {
//     uint64_t value;         // 64 bit value, send LSB first.
//     int n_block_radio;      // Size of transmitted block, all in bytes.
//     int k_data_radio;       // Size of transmitted data part.
//     int n_block_rs;         // Size of RS algorithm block.
//     int k_data_rs;          // Size of RS algorithm data part.
//     int itab;               // Index into Tab array.
// };
//
// static const struct correlation_tag_s tags[16] = {
//     /* Tag_00 */{ 0x566ED2717946107ELL,   0,   0,   0,   0, -1 },  //  Reserved
//
//     /* Tag_01 */{ 0xB74DB7DF8A532F3ELL, 255, 239, 255, 239, 0 },  //  RS(255, 239) 16-byte check value, 239 information bytes
//     /* Tag_02 */{ 0x26FF60A600CC8FDELL, 144, 128, 255, 239, 0 },  //  RS(144,128) - shortened RS(255, 239), 128 info bytes
//     /* Tag_03 */{ 0xC7DC0508F3D9B09ELL,  80,  64, 255, 239, 0 },  //  RS(80,64) - shortened RS(255, 239), 64 info bytes
//     /* Tag_04 */{ 0x8F056EB4369660EELL,  48,  32, 255, 239, 0 },  //  RS(48,32) - shortened RS(255, 239), 32 info bytes
//
//     /* Tag_05 */{ 0x6E260B1AC5835FAELL, 255, 223, 255, 223, 1 },  //  RS(255, 223) 32-byte check value, 223 information bytes
//     /* Tag_06 */{ 0xFF94DC634F1CFF4ELL, 160, 128, 255, 223, 1 },  //  RS(160,128) - shortened RS(255, 223), 128 info bytes
//     /* Tag_07 */{ 0x1EB7B9CDBC09C00ELL,  96,  64, 255, 223, 1 },  //  RS(96,64) - shortened RS(255, 223), 64 info bytes
//     /* Tag_08 */{ 0xDBF869BD2DBB1776LL,  64,  32, 255, 223, 1 },  //  RS(64,32) - shortened RS(255, 223), 32 info bytes
//
//     /* Tag_09 */{ 0x3ADB0C13DEAE2836LL, 255, 191, 255, 191, 2 },  //  RS(255, 191) 64-byte check value, 191 information bytes
//     /* Tag_0A */{ 0xAB69DB6A543188D6LL, 192, 128, 255, 191, 2 },  //  RS(192, 128) - shortened RS(255, 191), 128 info bytes
//     /* Tag_0B */{ 0x4A4ABEC4A724B796LL, 128,  64, 255, 191, 2 },  //  RS(128, 64) - shortened RS(255, 191), 64 info bytes
//
//     /* Tag_0C */{ 0x0293D578626B67E6LL,   0,   0,   0,   0, -1 },  //  Undefined
//     /* Tag_0D */{ 0xE3B0B0D6917E58A6LL,   0,   0,   0,   0, -1 },  //  Undefined
//     /* Tag_0E */{ 0x720267AF1BE1F846LL,   0,   0,   0,   0, -1 },  //  Undefined
//     /* Tag_0F */{ 0x93210201E8F4C706LL,   0,   0,   0,   0, -1 }   //  Undefined
// };

