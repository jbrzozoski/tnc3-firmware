#pragma once

// #include "HdlcFrame.hpp"
// #include "main.h"
// #include <cmsis_os.h>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <cassert>

//namespace mobilinkd { namespace tnc {

// TODO - Make these dynamic based on the parameter sets we are using
#define MAX_LOG_TABLE_SIZE 256
#define MAX_GENPOLY_SIZE  65

template <typename DTYPE>
struct rs_config {
    // DTYPE: data type being operated on (expected to be uint8_t)
    const uint16_t SYMSIZE; /* Bits per symbol (aka MM in some reference code) (expected to be 8) */
    const uint16_t GFPOLY; /* Field generator polynomial coefficients */
    const uint8_t FCR; /* First consecutive root, index form */
    const uint8_t PRIM; /* Primitive element, index form */
    const uint16_t NROOTS; /* Number of generator roots = number of parity symbols */
    const uint16_t NN = ((1<<SYMSIZE)-1); /* Symbols per block */
    const uint16_t A0 = NN;
    uint8_t alpha_to[MAX_LOG_TABLE_SIZE]; /* log lookup table */
    uint8_t index_of[MAX_LOG_TABLE_SIZE]; /* Antilog lookup table */
    uint8_t genpoly[MAX_GENPOLY_SIZE];    /* Generator polynomial */
    uint8_t iprim;       /* prim-th root of 1, index form */

    inline uint16_t modnn(uint16_t x) const {
        while (x >= NN) {
            x -= NN;
            x = (x >> SYMSIZE) + (x & NN);
        }
        return x;
    }

    rs_config(uint16_t symsize, uint16_t gfpoly, uint8_t fcr, uint8_t prim, uint16_t nroots)
    : SYMSIZE(symsize), GFPOLY(gfpoly), FCR(fcr), PRIM(prim), NROOTS(nroots), alpha_to(), index_of(), genpoly(), iprim(0)
    {
        int i, j, sr, root;

        assert(SYMSIZE <= (8 * sizeof(DTYPE)));
        assert(FCR < (1 << SYMSIZE));
        assert((PRIM != 0) && (PRIM < (1 << SYMSIZE)));
        assert(NROOTS < (1 << SYMSIZE));
        assert((sizeof(DTYPE) * (NN + 1)) <= MAX_LOG_TABLE_SIZE);
        assert((sizeof(DTYPE) * (NROOTS + 1)) <= MAX_GENPOLY_SIZE);

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

    //void encode(IoFrame* frame_input, IoFrame* fec_output) const
    void encode(DTYPE * input, DTYPE * output) const
    {
        int i, j;
        DTYPE feedback;
        DTYPE fec_work[NROOTS];

        // Input frame is expected to be of size NN-NROOTS
        // Generated output will be of size NROOTS

        // Clear out the FEC work area.
        // Every iteration of the loop we need to shift the full contents down by one
        // ex and most of the time also XOR every value.
        // The old design did the XOR first, then did a memmove.
        // This design does both the move and the XOR in one operation.
        memset(fec_work, 0, sizeof(fec_work));

        for(i = 0; i < (NN-NROOTS); i++) {
            feedback = index_of[input[i] ^ fec_work[0]];
            if (feedback != A0) {
                // Feedback is non-zero, so we XOR everything as we shift the contents
                for(j = 1; j < NROOTS; j++) {
                    fec_work[j-1] = fec_work[j] ^ alpha_to[modnn(feedback + genpoly[NROOTS-j])];
                }
                fec_work[NROOTS-1] = alpha_to[modnn(feedback + genpoly[0])];
            }
            else {
                // Feedback is zero, so we just shift the contents
                memmove(&fec_work[0], &fec_work[1], sizeof(DTYPE)*(NROOTS-1));
                fec_work[NROOTS-1] = 0;
            }
        }


        //fec_output.clear();
        //for (i = 0; i < NROOTS; i++) {
        //    bb.push_back(fec_work[i]);
        //}
        memcpy(output, &fec_work[0], sizeof(DTYPE)*NROOTS);
    }
};

struct fx25 {
    const unsigned FX25_BLOCK_SIZE = 255;

    enum RS_CONFIG_ENUM {
        RS_255_239 = 0,
        RS_255_223 = 1,
        RS_255_191 = 2,
        MAX_RS_CONFIG_ENUM
    };

    const struct rs_config<uint8_t> RS_CONFIGS[MAX_RS_CONFIG_ENUM] = {
        [RS_255_239] = rs_config(8, 0x11d, 1, 1, 16), // RS(255, 239) 16-byte check value, 239 information bytes
        [RS_255_223] = rs_config(8, 0x11d, 1, 1, 32), // RS(255, 223) 32-byte check value, 223 information bytes
        [RS_255_191] = rs_config(8, 0x11d, 1, 1, 64), // RS(255, 191) 64-byte check value, 191 information bytes
    };

    struct correlation_tag_s {
        uint64_t value;         // 64 bit value, send LSB first.
        int n_block_radio;      // Size of transmitted block, all in bytes.
        int k_data_radio;       // Size of transmitted data part.
        int n_block_rs;         // Size of RS algorithm block.
        int k_data_rs;          // Size of RS algorithm data part.
        enum RS_CONFIG_ENUM rs_config_index; // Index into RS_CONFIGS array.
    };

    const unsigned NUM_TAGS_DEFINED = 16;
    const unsigned CTAG_MIN = 0x01;
    const unsigned CTAG_MAX = 0x0B;
    const struct correlation_tag_s TAGS[NUM_TAGS_DEFINED] = {
        /* Tag_00 */{ 0x566ED2717946107ELL,   0,   0,   0,   0, -1 },  //  Reserved

        /* Tag_01 */{ 0xB74DB7DF8A532F3ELL, 255, 239, 255, 239, RS_255_239 },  //  RS(255, 239) 16-byte check value, 239 information bytes
        /* Tag_02 */{ 0x26FF60A600CC8FDELL, 144, 128, 255, 239, RS_255_239 },  //  RS(144,128) - shortened RS(255, 239), 128 info bytes
        /* Tag_03 */{ 0xC7DC0508F3D9B09ELL,  80,  64, 255, 239, RS_255_239 },  //  RS(80,64) - shortened RS(255, 239), 64 info bytes
        /* Tag_04 */{ 0x8F056EB4369660EELL,  48,  32, 255, 239, RS_255_239 },  //  RS(48,32) - shortened RS(255, 239), 32 info bytes

        /* Tag_05 */{ 0x6E260B1AC5835FAELL, 255, 223, 255, 223, RS_255_223 },  //  RS(255, 223) 32-byte check value, 223 information bytes
        /* Tag_06 */{ 0xFF94DC634F1CFF4ELL, 160, 128, 255, 223, RS_255_223 },  //  RS(160,128) - shortened RS(255, 223), 128 info bytes
        /* Tag_07 */{ 0x1EB7B9CDBC09C00ELL,  96,  64, 255, 223, RS_255_223 },  //  RS(96,64) - shortened RS(255, 223), 64 info bytes
        /* Tag_08 */{ 0xDBF869BD2DBB1776LL,  64,  32, 255, 223, RS_255_223 },  //  RS(64,32) - shortened RS(255, 223), 32 info bytes

        /* Tag_09 */{ 0x3ADB0C13DEAE2836LL, 255, 191, 255, 191, RS_255_191 },  //  RS(255, 191) 64-byte check value, 191 information bytes
        /* Tag_0A */{ 0xAB69DB6A543188D6LL, 192, 128, 255, 191, RS_255_191 },  //  RS(192, 128) - shortened RS(255, 191), 128 info bytes
        /* Tag_0B */{ 0x4A4ABEC4A724B796LL, 128,  64, 255, 191, RS_255_191 },  //  RS(128, 64) - shortened RS(255, 191), 64 info bytes

        /* Tag_0C */{ 0x0293D578626B67E6LL,   0,   0,   0,   0, -1 },  //  Undefined
        /* Tag_0D */{ 0xE3B0B0D6917E58A6LL,   0,   0,   0,   0, -1 },  //  Undefined
        /* Tag_0E */{ 0x720267AF1BE1F846LL,   0,   0,   0,   0, -1 },  //  Undefined
        /* Tag_0F */{ 0x93210201E8F4C706LL,   0,   0,   0,   0, -1 }   //  Undefined
    }

    fx25() {
        // Verify integrity of tables and assumptions.
        // This also does a quick check for the popcount function.
        for (int j = 0; j < NUM_TAGS_DEFINED; j++) {
            for (int k = 0; k < NUM_TAGS_DEFINED; k++) {
                if (j == k) {
                    assert(__builtin_popcountll(TAGS[j].value ^ TAGS[k].value) == 0);
                } else {
                    assert(__builtin_popcountll(TAGS[j].value ^ TAGS[k].value) == 32);
                }
            }
        }

        for (int j = CTAG_MIN; j <= CTAG_MAX; j++) {
            assert(TAGS[j].n_block_radio - TAGS[j].k_data_radio == RS_CONFIGS[TAGS[j].rs_config_index].NROOTS);
            assert(TAGS[j].n_block_rs - TAGS[j].k_data_rs == RS_CONFIGS[TAGS[j].rs_config_index].NROOTS);
            assert(TAGS[j].n_block_rs == FX25_BLOCK_SIZE);
        }

        assert(pick_mode(1, 239) == 1);
        assert(pick_mode(1, 240) == -1);

        assert(pick_mode(5, 223) == 5);
        assert(pick_mode(5, 224) == -1);

        assert(pick_mode(9, 191) == 9);
        assert(pick_mode(9, 192) == -1);

        assert(pick_mode(16, 32) == 4);
        assert(pick_mode(16, 64) == 3);
        assert(pick_mode(16, 128) == 2);
        assert(pick_mode(16, 239) == 1);
        assert(pick_mode(16, 240) == -1);

        assert(pick_mode(32, 32) == 8);
        assert(pick_mode(32, 64) == 7);
        assert(pick_mode(32, 128) == 6);
        assert(pick_mode(32, 223) == 5);
        assert(pick_mode(32, 234) == -1);

        assert(pick_mode(64, 64) == 11);
        assert(pick_mode(64, 128) == 10);
        assert(pick_mode(64, 191) == 9);
        assert(pick_mode(64, 192) == -1);
    }

    int pick_mode(int fx_mode, int dlen) const {
        // Does it look like fx_mode is a direct index to a usable TAGS entry?
        if (fx_mode >= CTAG_MIN && fx_mode <= CTAG_MAX) {
            // If so, just return the index back or fail based on whether the dlen will fit
            return ((dlen <= TAGS[fx_mode].k_data_radio) ? fx_mode : -1);
        }

        // Otherwise, assume fx_mode is a request for a number of NROOTS,
        // and find the smallest mode with that number that will fit dlen data
        for (int k = CTAG_MAX; k >= CTAG_MIN; k--) {
            if ((fx_mode == RS_CONFIGS[TAGS[k].rs_config_index].NROOTS) &&
                (dlen <= TAGS[fx_mode].k_data_radio)) {
                return (k);
            }
        }
        return (-1);
    }

    // void encode(IoFrame *frame) const {
    // }
}

