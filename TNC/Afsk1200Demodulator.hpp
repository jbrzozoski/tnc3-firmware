// Copyright 2020 Rob Riggs <rob@mobilinkd.com>
// All rights reserved.

#pragma once

#include "Demodulator.hpp"
#include "AfskDemodulator.hpp"
#include "FirFilter.hpp"
#include "FilterCoefficients.hpp"
#include "KissHardware.hpp"
#include "HdlcFrame.hpp"

namespace mobilinkd { namespace tnc {

struct Afsk1200Demodulator : IDemodulator
{
    /*
     * Generated with Scipy Filter, 152 coefficients, 1100-2300Hz bandpass,
     * Hann window, starting and ending 0 value coefficients removed.
     *
     * np.array(
     *  firwin2(152,
     *      [
     *          0.0,
     *          1000.0/(sample_rate/2),
     *          1100.0/(sample_rate/2),
     *          2350.0/(sample_rate/2),
     *          2500.0/(sample_rate/2),
     *          1.0
     *      ],
     *      [0,0,1,1,0,0],
     *      antisymmetric = False,
     *      window='hann') * 32768,
     *  dtype=int)[10:-10]
     */
    static constexpr size_t FILTER_TAP_NUM = 132;
    static constexpr uint32_t ADC_BLOCK_SIZE = afsk1200::ADC_BUFFER_SIZE;
    static constexpr uint32_t SAMPLE_RATE = 26400;
    static_assert(audio::ADC_BUFFER_SIZE >= ADC_BLOCK_SIZE);

    using audio_filter_t = Q15FirFilter<ADC_BLOCK_SIZE, FILTER_TAP_NUM>;

    static const q15_t bpf_coeffs[FILTER_TAP_NUM];


    static afsk1200::emphasis_filter_type filter_1;
    static afsk1200::emphasis_filter_type filter_2;
    static afsk1200::emphasis_filter_type filter_3;

    static afsk1200::Demodulator demod1;
    static afsk1200::Demodulator demod2;
    static afsk1200::Demodulator demod3;

    audio_filter_t demod_filter;
    uint16_t last_fcs{0};
    uint32_t last_counter{0};
    uint32_t counter{0};
    bool locked_{false};

    virtual ~Afsk1200Demodulator() {}

    void start() override
    {
        INFO("Setting 48MHz SysClock.");
        SysClock48();

        // rx_twist is 6dB for discriminator input and 0db for de-emphasized input.
        auto twist = kiss::settings().rx_twist;

        filter_1.init(*filter::fir::AfskFilters[twist + 3]);
        filter_2.init(*filter::fir::AfskFilters[twist + 6]);
        filter_3.init(*filter::fir::AfskFilters[twist + 9]);

        last_fcs = 0;
        last_counter = 0;
        counter = 0;

        demod_filter.init(bpf_coeffs);

        hadc1.Init.OversamplingMode = ENABLE;
        if (HAL_ADC_Init(&hadc1) != HAL_OK)
        {
            CxxErrorHandler();
        }

        ADC_ChannelConfTypeDef sConfig;

        sConfig.Channel = AUDIO_IN;
        sConfig.Rank = ADC_REGULAR_RANK_1;
        sConfig.SingleDiff = ADC_SINGLE_ENDED;
        sConfig.SamplingTime = ADC_SAMPLETIME_12CYCLES_5;
        sConfig.OffsetNumber = ADC_OFFSET_NONE;
        sConfig.Offset = 0;
        if (HAL_ADC_ConfigChannel(&hadc1, &sConfig) != HAL_OK)
            CxxErrorHandler();

        startADC(1817, ADC_BLOCK_SIZE);
    }

    void stop() override
    {
        stopADC();
        locked_ = false;
    }

    hdlc::IoFrame* operator()(const q15_t* samples) override
    {
        hdlc::IoFrame* result = nullptr;

        q15_t* filtered = demod_filter(const_cast<q15_t* >(samples));

        ++counter;

#if 1
        auto frame1 = demod1(filtered, ADC_BLOCK_SIZE);
        if (frame1)
        {
            if (frame1->fcs() != last_fcs or counter > last_counter + 2)
            {
                last_fcs = frame1->fcs();
                last_counter = counter;
                result = frame1;
            }
            else
            {
                hdlc::release (frame1);
            }
        }
#endif

#if 1
        auto frame2 = demod2(filtered, ADC_BLOCK_SIZE);
        if (frame2)
        {
            if (frame2->fcs() != last_fcs or counter > last_counter + 2)
            {
                last_fcs = frame2->fcs();
                last_counter = counter;
                result = frame2;
            }
            else
            {
                hdlc::release(frame2);
            }
        }
#endif

#if 1
        auto frame3 = demod3(filtered, ADC_BLOCK_SIZE);
        if (frame3)
        {
            if (frame3->fcs() != last_fcs or counter > last_counter + 2)
            {
                last_fcs = frame3->fcs();
                last_counter = counter;
                result = frame3;
            }
            else
            {
                hdlc::release(frame3);
            }
        }
#endif
        locked_ = demod1.locked() or demod2.locked() or demod3.locked();
        return result;
    }

    float readTwist() override;

    bool locked() const override
    {
        return locked_;
    }

    size_t size() const override
    {
        return ADC_BLOCK_SIZE;
    }
};

}} // mobilinkd::tnc
