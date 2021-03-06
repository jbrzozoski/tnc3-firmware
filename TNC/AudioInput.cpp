// Copyright 2018-2019 Rob Riggs <rob@mobilinkd.com>
// All rights reserved.

#include "AudioInput.hpp"
#include "AfskDemodulator.hpp"
#include "AudioLevel.hpp"
#include "Log.h"
#include "KissHardware.hpp"
#include "GPIO.hpp"
#include "HdlcFrame.hpp"
#include "memory.hpp"
#include "IirFilter.hpp"
#include "FilterCoefficients.hpp"
#include "PortInterface.hpp"
#include "Goertzel.h"
#include "DCD.h"

#include "arm_math.h"
#include "stm32l4xx_hal.h"

#include <algorithm>
#include <numeric>
#include <cstring>
#include <cstdint>
#include <atomic>

extern osMessageQId ioEventQueueHandle;

extern "C" void SystemClock_Config(void);

// 1kB
typedef mobilinkd::tnc::memory::Pool<
    8, mobilinkd::tnc::audio::ADC_BUFFER_SIZE * 2> adc_pool_type;
adc_pool_type adcPool;

// DMA Conversion first half complete.
extern "C" void HAL_ADC_ConvHalfCpltCallback(ADC_HandleTypeDef*) {
    using namespace mobilinkd::tnc::audio;

    auto block = adcPool.allocate();
    if (!block) return;
    memmove(block->buffer, adc_buffer, ADC_BUFFER_SIZE * 2);
    auto status = osMessagePut(adcInputQueueHandle, (uint32_t) block, 0);
    if (status != osOK) adcPool.deallocate(block);
}

// DMA Conversion second half complete.
extern "C" void HAL_ADC_ConvCpltCallback(ADC_HandleTypeDef*) {
    using namespace mobilinkd::tnc::audio;

    auto block = adcPool.allocate();
    if (!block) return;
    memmove(block->buffer, adc_buffer + DMA_TRANSFER_SIZE, ADC_BUFFER_SIZE * 2);
    auto status = osMessagePut(adcInputQueueHandle, (uint32_t) block, 0);
    if (status != osOK) adcPool.deallocate(block);
}

extern "C" void HAL_ADC_ErrorCallback(ADC_HandleTypeDef* /* hadc */) {
    using namespace mobilinkd::tnc::audio;

    // __HAL_ADC_CLEAR_FLAG(hadc, (ADC_FLAG_EOC | ADC_FLAG_EOS | ADC_FLAG_OVR));
    // HAL_DMA_Start(hadc->DMA_Handle, (uint32_t)&hadc->Instance->DR, (uint32_t)adc_buffer, ADC_BUFFER_SIZE * 2);
}

extern "C" void startAudioInputTask(void const*) {

    using namespace mobilinkd::tnc::audio;
    DEBUG("startAudioInputTask");

    adcPool.init();

    uint8_t adcState = mobilinkd::tnc::audio::IDLE;

    while (true) {
        osEvent event = osMessageGet(audioInputQueueHandle, osWaitForever);
        if (event.status != osEventMessage) continue;
        adcState = event.value.v;

        switch (adcState) {
        case STOPPED:
            DEBUG("STOPPED");
            // stop();
            break;
        case DEMODULATOR:
            DEBUG("DEMODULATOR");
            demodulatorTask();
            break;
        case STREAM_AMPLIFIED_INPUT_LEVEL:
            DEBUG("STREAM_AMPLIFIED_INPUT_LEVEL");
            streamAmplifiedInputLevels();
            break;
        case POLL_AMPLIFIED_INPUT_LEVEL:
            DEBUG("POLL_AMPLIFIED_INPUT_LEVEL");
            pollAmplifiedInputLevel();
            break;
        case POLL_BATTERY_LEVEL:
            DEBUG("POLL_BATTERY_LEVEL");
            pollBatteryLevel();
            break;
        case POLL_TWIST_LEVEL:
            DEBUG("POLL_TWIST_LEVEL");
            pollInputTwist();
            break;
        case STREAM_AVERAGE_TWIST_LEVEL:
            DEBUG("STREAM_AVERAGE_TWIST_LEVEL");
            // streamAverageInputTwist();
            break;
        case STREAM_INSTANT_TWIST_LEVEL:
            DEBUG("STREAM_INSTANT_TWIST_LEVEL");
            // streamInstantInputTwist();
            break;
        case AUTO_ADJUST_INPUT_LEVEL:
            DEBUG("AUTO_ADJUST_INPUT_LEVEL");
            autoAudioInputLevel();
            break;
        case CONFIGURE_INPUT_LEVELS:
            DEBUG("CONFIGURE_INPUT_LEVELS");
            setAudioInputLevels();
            break;
        case UPDATE_SETTINGS:
            DEBUG("UPDATE_SETTINGS");
            setAudioInputLevels();
            break;
        case IDLE:
            DEBUG("IDLE");
            break;
        default:
            break;
        }
    }
}

namespace mobilinkd { namespace tnc { namespace audio {

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
constexpr size_t FILTER_TAP_NUM = 132;
const q15_t bpf_coeffs[] = {
    4,     0,    -5,   -10,   -13,   -12,    -9,    -4,    -2,    -4,   -12,   -26,
  -41,   -52,   -51,   -35,    -3,    39,    83,   117,   131,   118,    83,    36,
   -6,   -32,   -30,    -3,    36,    67,    66,    19,   -74,  -199,  -323,  -408,
 -421,  -344,  -187,    17,   218,   364,   417,   369,   247,   106,    14,    26,
  166,   407,   676,   865,   866,   605,    68,  -675, -1484, -2171, -2547, -2471,
-1895,  -882,   394,  1692,  2747,  3337,  3337,  2747,  1692,   394,  -882, -1895,
-2471, -2547, -2171, -1484,  -675,    68,   605,   866,   865,   676,   407,   166,
   26,    14,   106,   247,   369,   417,   364,   218,    17,  -187,  -344,  -421,
 -408,  -323,  -199,   -74,    19,    66,    67,    36,    -3,   -30,   -32,    -6,
   36,    83,   118,   131,   117,    83,    39,    -3,   -35,   -51,   -52,   -41,
  -26,   -12,    -4,    -2,    -4,    -9,   -12,   -13,   -10,    -5,     0,     4,
};

uint32_t adc_buffer[ADC_BUFFER_SIZE];       // Two samples per element.

typedef Q15FirFilter<ADC_BUFFER_SIZE, FILTER_TAP_NUM> audio_filter_type;

audio_filter_type audio_filter;

mobilinkd::tnc::afsk1200::emphasis_filter_type filter_1;
mobilinkd::tnc::afsk1200::emphasis_filter_type filter_2;
mobilinkd::tnc::afsk1200::emphasis_filter_type filter_3;

mobilinkd::tnc::afsk1200::Demodulator& getDemod1(const TFirCoefficients<9>& f) __attribute__((noinline));
mobilinkd::tnc::afsk1200::Demodulator& getDemod2(const TFirCoefficients<9>& f) __attribute__((noinline));
mobilinkd::tnc::afsk1200::Demodulator& getDemod3(const TFirCoefficients<9>& f) __attribute__((noinline));

mobilinkd::tnc::afsk1200::Demodulator& getDemod1(const TFirCoefficients<9>& f) {
    filter_1.init(f);
    static mobilinkd::tnc::afsk1200::Demodulator instance(26400, filter_1);
    return instance;
}

mobilinkd::tnc::afsk1200::Demodulator& getDemod2(const TFirCoefficients<9>& f) {
    filter_2.init(f);
    static mobilinkd::tnc::afsk1200::Demodulator instance(26400, filter_2);
    return instance;
}

mobilinkd::tnc::afsk1200::Demodulator& getDemod3(const TFirCoefficients<9>& f) {
    filter_3.init(f);
    static mobilinkd::tnc::afsk1200::Demodulator instance(26400, filter_3);
    return instance;
}

q15_t normalized[ADC_BUFFER_SIZE];

void demodulatorTask() {

    DEBUG("enter demodulatorTask");

    audio_filter.init(bpf_coeffs);

    // rx_twist is 6dB for discriminator input and 0db for de-emphasized input.
    auto twist = kiss::settings().rx_twist;

    mobilinkd::tnc::afsk1200::Demodulator& demod1 = getDemod1(*filter::fir::AfskFilters[twist + 3]);
    mobilinkd::tnc::afsk1200::Demodulator& demod2 = getDemod2(*filter::fir::AfskFilters[twist + 6]);
    mobilinkd::tnc::afsk1200::Demodulator& demod3 = getDemod3(*filter::fir::AfskFilters[twist + 9]);

    startADC(AUDIO_IN);

    mobilinkd::tnc::hdlc::IoFrame* frame = 0;

    uint16_t last_fcs = 0;
    uint32_t last_counter = 0;
    uint32_t counter = 0;

    bool dcd_status{false};

    while (true) {
        osEvent peek = osMessagePeek(audioInputQueueHandle, 0);
        if (peek.status == osEventMessage) break;

        osEvent evt = osMessageGet(adcInputQueueHandle, osWaitForever);
        if (evt.status != osEventMessage) {
            continue;
        }

        ++counter;

        auto block = (adc_pool_type::chunk_type*) evt.value.p;
        auto samples = (int16_t*) block->buffer;

        arm_offset_q15(samples, 0 - virtual_ground, normalized, ADC_BUFFER_SIZE);
        adcPool.deallocate(block);
        q15_t* audio = audio_filter(normalized);

#if 1
        frame = demod1(audio, ADC_BUFFER_SIZE);
        if (frame) {
            if (frame->fcs() != last_fcs or counter > last_counter + 2) {
                auto save_fcs = frame->fcs();
                if (osMessagePut(ioEventQueueHandle, (uint32_t) frame, 1) == osOK) {
                  last_fcs = save_fcs;
                  last_counter = counter;
                } else {
                  hdlc::release(frame);
                }
            }
            else {
                hdlc::release(frame);
            }
        }
#endif

#if 1
        frame = demod2(audio, ADC_BUFFER_SIZE);
        if (frame) {
          if (frame->fcs() != last_fcs or counter > last_counter + 2) {
              auto save_fcs = frame->fcs();
              if (osMessagePut(ioEventQueueHandle, (uint32_t) frame, 1) == osOK) {
                last_fcs = save_fcs;
                last_counter = counter;
              } else {
                hdlc::release(frame);
              }
          }
          else {
              hdlc::release(frame);
          }
        }
#endif

#if 1
        frame = demod3(audio, ADC_BUFFER_SIZE);
        if (frame) {
          if (frame->fcs() != last_fcs or counter > last_counter + 2) {
              auto save_fcs = frame->fcs();
              if (osMessagePut(ioEventQueueHandle, (uint32_t) frame, 1) == osOK) {
                last_fcs = save_fcs;
                last_counter = counter;
              } else {
                hdlc::release(frame);
              }
          }
          else {
              hdlc::release(frame);
          }
        }
#endif
        bool new_dcd_status = demod1.locked() or demod2.locked() or demod3.locked();
        if (new_dcd_status xor dcd_status) {
            dcd_status = new_dcd_status;
            if (dcd_status) {
                dcd_on();
            } else {
                dcd_off();
            }
        }
    }

    stopADC();
    dcd_off();
    DEBUG("exit demodulatorTask");
}


void streamLevels(uint32_t channel, uint8_t cmd) {

    // Stream out Vpp, Vavg, Vmin, Vmax as four 16-bit values, left justified.

    uint8_t data[9];
    INFO("streamLevels: start");
    startADC(channel);

    while (true) {
        osEvent peek = osMessagePeek(audioInputQueueHandle, 0);
        if (peek.status == osEventMessage) break;

        uint16_t count = 0;
        uint32_t accum = 0;
        uint16_t vmin = std::numeric_limits<uint16_t>::max();
        uint16_t vmax = std::numeric_limits<uint16_t>::min();

        while (count < 2640) {
            osEvent evt = osMessageGet(adcInputQueueHandle, osWaitForever);
            if (evt.status != osEventMessage) continue;

            count += ADC_BUFFER_SIZE;

            auto block = (adc_pool_type::chunk_type*) evt.value.p;
            auto start =  (uint16_t*) block->buffer;
            auto end = start + ADC_BUFFER_SIZE;

            vmin = std::min(vmin, *std::min_element(start, end));
            vmax = std::max(vmax, *std::max_element(start, end));
            accum = std::accumulate(start, end, accum);

            adcPool.deallocate(block);
        }

        uint16_t pp = (vmax - vmin) << 4;
        uint16_t avg = (accum / count) << 4;
        vmin <<= 4;
        vmax <<= 4;

        data[0] = cmd;
        data[1] = (pp >> 8) & 0xFF;   // Vpp
        data[2] = (pp & 0xFF);
        data[3] = (avg >> 8) & 0xFF;  // Vavg (DC level)
        data[4] = (avg & 0xFF);
        data[5] = (vmin >> 8) & 0xFF;  // Vmin
        data[6] = (vmin & 0xFF);
        data[7] = (vmax >> 8) & 0xFF;  // Vmax
        data[8] = (vmax & 0xFF);

        ioport->write(data, 9, 6, 10);
    }

    stopADC();
    DEBUG("exit streamLevels");
}

levels_type readLevels(uint32_t channel, uint32_t samples) {

    DEBUG("enter readLevels");

    // Return Vpp, Vavg, Vmin, Vmax as four 16-bit values, right justified.

    uint16_t count = 0;
    uint32_t accum = 0;
    uint16_t vmin = std::numeric_limits<uint16_t>::max();
    uint16_t vmax = std::numeric_limits<uint16_t>::min();

    INFO("readLevels: start");
    startADC(channel);

    while (count < samples) {

        osEvent evt = osMessageGet(adcInputQueueHandle, osWaitForever);
        if (evt.status != osEventMessage) continue;

        auto block = (adc_pool_type::chunk_type*) evt.value.p;
        auto start =  (uint16_t*) block->buffer;
        auto end = start + ADC_BUFFER_SIZE;

        vmin = std::min(vmin, *std::min_element(start, end));
        vmax = std::max(vmax, *std::max_element(start, end));
        accum = std::accumulate(start, end, accum);

        adcPool.deallocate(block);

        count += ADC_BUFFER_SIZE;
    }

    stopADC();

    uint16_t pp = vmax - vmin;
    uint16_t avg = accum / count;
    DEBUG("exit readLevels");

    return levels_type(pp, avg, vmin, vmax);
}


/**
 * This provides 100Hz resolution to the Goerztel filter.
 */
constexpr uint32_t TWIST_SAMPLE_SIZE = 88;

/*
 * Return twist as a the difference in dB between mark and space.  The
 * expected values are about 0dB for discriminator output and about 5.5dB
 * for de-emphasized audio.
 */
float readTwist()
{

  DEBUG("enter readTwist");
  constexpr uint32_t channel = AUDIO_IN;

  float g1200 = 0.0f;
  float g2200 = 0.0f;

  GoertzelFilter<TWIST_SAMPLE_SIZE, SAMPLE_RATE> gf1200(1200.0, 0);
  GoertzelFilter<TWIST_SAMPLE_SIZE, SAMPLE_RATE> gf2200(2200.0, 0);

  const uint32_t AVG_SAMPLES = 20;

  startADC(channel);

  for (uint32_t i = 0; i != AVG_SAMPLES; ++i)
  {
    uint32_t count = 0;
    while (count < TWIST_SAMPLE_SIZE)
    {
        osEvent evt = osMessageGet(adcInputQueueHandle, osWaitForever);
        if (evt.status != osEventMessage) continue;

        auto block = (adc_pool_type::chunk_type*) evt.value.p;
        uint16_t* data =  (uint16_t*) block->buffer;
        gf1200(data, ADC_BUFFER_SIZE);
        gf2200(data, ADC_BUFFER_SIZE);

        adcPool.deallocate(block);

        count += ADC_BUFFER_SIZE;
    }

    g1200 += (gf1200 / count);
    g2200 += (gf2200 / count);

    gf1200.reset();
    gf2200.reset();
  }

  stopADC();

  g1200 = 10.0f * log10f(g1200 / AVG_SAMPLES);
  g2200 = 10.0f * log10f(g2200 / AVG_SAMPLES);

  auto result = g1200 - g2200;

  INFO("Twist = %d / 100 (%d - %d)", int(result*100), int(g1200), int(g2200));

  DEBUG("exit readTwist");
  return result;
}

/*
 * Get the input twist level as a pair of numbers -- the relative dB
 * level of the Bell 202 mark and space tones.
 *
 * This is intended to measure noise levels on an empty channel.
 *
 * When de-emphasis is applied, the noise at 1200Hz will be about 5.5dB
 * higher than at 2200Hz.  When de-emphasis is not applied (discriminator
 * output), the levels should be about the same.
 *
 * This is used to adjust the demodulator filters so that the proper
 * input twist is applied to the signal.  In general, properly modulated
 * signals are expected to be pre-emphasized so that they are equal
 * when de-emphasis is applied.
 *
 * If no de-emphasis is detected, the de-emphasis has to be applied in
 * the demodulator.
 *
 * This takes about 5 seconds to complete as it averages 100 50ms samples
 * to get a reasonable sampling of the noise.
 */
void pollInputTwist()
{
    DEBUG("enter pollInputTwist");
    constexpr uint32_t channel = AUDIO_IN;

    float g1200 = 0.0f;
    float g2200 = 0.0f;

    GoertzelFilter<TWIST_SAMPLE_SIZE, SAMPLE_RATE> gf1200(1200.0);
    GoertzelFilter<TWIST_SAMPLE_SIZE, SAMPLE_RATE> gf2200(2200.0);

    const uint32_t AVG_SAMPLES = 100;

    startADC(channel);

    for (uint32_t i = 0; i != AVG_SAMPLES; ++i) {

      uint32_t count = 0;
      while (count < TWIST_SAMPLE_SIZE) {

          osEvent evt = osMessageGet(adcInputQueueHandle, osWaitForever);
          if (evt.status != osEventMessage) continue;

          count += ADC_BUFFER_SIZE;

          auto block = (adc_pool_type::chunk_type*) evt.value.p;
          uint16_t* data =  (uint16_t*) block->buffer;
          gf1200(data, ADC_BUFFER_SIZE);
          gf2200(data, ADC_BUFFER_SIZE);

          adcPool.deallocate(block);
      }

      g1200 += 10.0 * log10(gf1200);
      g2200 += 10.0 * log10(gf2200);

      gf1200.reset();
      gf2200.reset();
    }

    stopADC();

    DEBUG("pollInputTwist: MARK=%d, SPACE=%d (x100)",
      int(g1200 * 100.0 / AVG_SAMPLES), int(g2200 * 100.0 / AVG_SAMPLES));

    int16_t g1200i = int16_t(g1200 * 256 / AVG_SAMPLES);
    int16_t g2200i = int16_t(g2200 * 256 / AVG_SAMPLES);

    uint8_t buffer[5];
    buffer[0] = kiss::hardware::POLL_INPUT_TWIST;
    buffer[1] = (g1200i >> 8) & 0xFF;
    buffer[2] = g1200i & 0xFF;
    buffer[3] = (g2200i >> 8) & 0xFF;
    buffer[4] = g2200i & 0xFF;

    ioport->write(buffer, 5, 6, 10);

    DEBUG("exit pollInputTwist");
}

#if 0
void streamAverageInputTwist()
{
    DEBUG("enter streamAverageInputTwist");

    constexpr uint32_t channel = AUDIO_IN;

    startADC(channel);

    uint32_t acount = 0;
    float g700 = 0.0f;
    float g1200 = 0.0f;
    float g1700 = 0.0f;
    float g2200 = 0.0f;
    float g2700 = 0.0f;

    GoertzelFilter<TWIST_SAMPLE_SIZE, 26400> gf700(700.0);
    GoertzelFilter<TWIST_SAMPLE_SIZE, 26400> gf1200(1200.0);
    GoertzelFilter<TWIST_SAMPLE_SIZE, 26400> gf1700(1700.0);
    GoertzelFilter<TWIST_SAMPLE_SIZE, 26400> gf2200(2200.0);
    GoertzelFilter<TWIST_SAMPLE_SIZE, 26400> gf2700(2700.0);

    while (true) {
      osEvent peek = osMessagePeek(audioInputQueueHandle, 0);
      if (peek.status == osEventMessage) break;

      acount++;
      uint32_t count = 0;
      while (count < TWIST_SAMPLE_SIZE) {

          osEvent evt = osMessageGet(adcInputQueueHandle, osWaitForever);
          if (evt.status != osEventMessage) continue;

          count += ADC_BUFFER_SIZE;

          auto block = (adc_pool_type::chunk_type*) evt.value.v;
          uint16_t* data =  (uint16_t*) block->buffer;
          gf700(data, ADC_BUFFER_SIZE);
          gf1200(data, ADC_BUFFER_SIZE);
          gf1700(data, ADC_BUFFER_SIZE);
          gf2200(data, ADC_BUFFER_SIZE);
          gf2700(data, ADC_BUFFER_SIZE);

          adcPool.deallocate(block);
      }

      g700 += 10.0 * log10(gf700);
      g1200 += 10.0 * log10(gf1200);
      g1700 += 10.0 * log10(gf1700);
      g2200 += 10.0 * log10(gf2200);
      g2700 += 10.0 * log10(gf2700);

      char* buffer = 0;
      // @TODO: Make re-entrant printf work (or convert to fixed-point).
      int len = asiprintf_r(
        &buffer,
        "_%f, %f, %f, %f, %f\r\n",
        g700 / acount,
        g1200 / acount,
        g1700 / acount,
        g2200 / acount,
        g2700 / acount);

      if (len > 0) {
        buffer[0] = kiss::hardware::POLL_INPUT_TWIST;
        ioport->write((uint8_t*)buffer, len - 1, 6, 10);
        free(buffer);
      }

      gf700.reset();
      gf1200.reset();
      gf1700.reset();
      gf2200.reset();
      gf2700.reset();
    }

    stopADC();
    DEBUG("exit streamAverageInputTwist");
}

void streamInstantInputTwist()
{
    DEBUG("enter streamInstantInputTwist");

    constexpr uint32_t channel = AUDIO_IN;

    startADC(channel);

    GoertzelFilter<TWIST_SAMPLE_SIZE, 26400> gf700(700.0);
    GoertzelFilter<TWIST_SAMPLE_SIZE, 26400> gf1200(1200.0);
    GoertzelFilter<TWIST_SAMPLE_SIZE, 26400> gf1700(1700.0);
    GoertzelFilter<TWIST_SAMPLE_SIZE, 26400> gf2200(2200.0);
    GoertzelFilter<TWIST_SAMPLE_SIZE, 26400> gf2700(2700.0);

    while (true) {
      osEvent peek = osMessagePeek(audioInputQueueHandle, 0);
      if (peek.status == osEventMessage) break;

      uint32_t count = 0;
      while (count < TWIST_SAMPLE_SIZE) {

          osEvent evt = osMessageGet(adcInputQueueHandle, osWaitForever);
          if (evt.status != osEventMessage) continue;

          count += ADC_BUFFER_SIZE;

          auto block = (adc_pool_type::chunk_type*) evt.value.v;
          uint16_t* data =  (uint16_t*) block->buffer;
          gf700(data, ADC_BUFFER_SIZE);
          gf1200(data, ADC_BUFFER_SIZE);
          gf1700(data, ADC_BUFFER_SIZE);
          gf2200(data, ADC_BUFFER_SIZE);
          gf2700(data, ADC_BUFFER_SIZE);

          adcPool.deallocate(block);
      }

      char* buffer = 0;
      // @TODO: Make re-entrant printf work (or convert to fixed-point).
      int len = asiprintf_r(
        &buffer,
        "_%f, %f, %f, %f, %f\r\n",
        10.0 * log10(gf700),
        10.0 * log10(gf1200),
        10.0 * log10(gf1700),
        10.0 * log10(gf2200),
        10.0 * log10(gf2700));

      if (len > 0) {
        buffer[0] = kiss::hardware::POLL_INPUT_TWIST;
        ioport->write((uint8_t*)buffer, len - 1, 6, 10);
        free(buffer);
      }

      gf700.reset();
      gf1200.reset();
      gf1700.reset();
      gf2200.reset();
      gf2700.reset();
    }

    stopADC();
    DEBUG("exit streamInstantInputTwist");
}
#endif

void streamAmplifiedInputLevels() {
    DEBUG("enter streamAmplifiedInputLevels");
    streamLevels(AUDIO_IN, kiss::hardware::POLL_INPUT_LEVEL);
    DEBUG("exit streamAmplifiedInputLevels");
}

void pollAmplifiedInputLevel() {
    DEBUG("enter pollAmplifiedInputLevel");

    uint16_t Vpp, Vavg, Vmin, Vmax;
    std::tie(Vpp, Vavg, Vmin, Vmax) = readLevels(AUDIO_IN);

    Vpp <<= 4;
    Vavg <<= 4;
    Vmin <<= 4;
    Vmax <<= 4;

    uint8_t data[9];
    data[0] = kiss::hardware::POLL_INPUT_LEVEL;
    data[1] = (Vpp >> 8) & 0xFF;   // Vpp
    data[2] = (Vpp & 0xFF);
    data[3] = (Vavg >> 8) & 0xFF;  // Vavg (DC level)
    data[4] = (Vavg & 0xFF);
    data[5] = (Vmin >> 8) & 0xFF;  // Vmin
    data[6] = (Vmin & 0xFF);
    data[7] = (Vmax >> 8) & 0xFF;  // Vmax
    data[8] = (Vmax & 0xFF);

    ioport->write(data, 9, 6, 10);
    DEBUG("exit pollAmplifiedInputLevel");
}

void pollBatteryLevel() {
    DEBUG("enter pollBatteryLevel");

    ADC_ChannelConfTypeDef sConfig;

    sConfig.Channel = ADC_CHANNEL_VREFINT;
    sConfig.Rank = ADC_REGULAR_RANK_1;
    sConfig.SingleDiff = ADC_SINGLE_ENDED;
    sConfig.SamplingTime = ADC_SAMPLETIME_247CYCLES_5;
    sConfig.OffsetNumber = ADC_OFFSET_NONE;
    sConfig.Offset = 0;
    if (HAL_ADC_ConfigChannel(&hadc1, &sConfig) != HAL_OK)
        CxxErrorHandler();

    htim6.Init.Period = 48000;
    if (HAL_TIM_Base_Init(&htim6) != HAL_OK) CxxErrorHandler();

    if (HAL_TIM_Base_Start(&htim6) != HAL_OK)
        CxxErrorHandler();

    if (HAL_ADC_Start(&hadc1) != HAL_OK) CxxErrorHandler();
    if (HAL_ADC_PollForConversion(&hadc1, 3) != HAL_OK) CxxErrorHandler();
    auto vrefint = HAL_ADC_GetValue(&hadc1);
    if (HAL_ADC_Stop(&hadc1) != HAL_OK) CxxErrorHandler();

    // Disable battery charging while measuring battery voltage.
    auto usb_ce = gpio::USB_CE::get();
    gpio::USB_CE::on();

    gpio::BAT_DIVIDER::off();
    HAL_Delay(1);

    sConfig.Channel = ADC_CHANNEL_15;
    if (HAL_ADC_ConfigChannel(&hadc1, &sConfig) != HAL_OK)
        CxxErrorHandler();

    uint32_t vbat = 0;
    if (HAL_ADC_Start(&hadc1) != HAL_OK) CxxErrorHandler();
    for (size_t i = 0; i != 8; ++i)
    {
        if (HAL_ADC_PollForConversion(&hadc1, 1) != HAL_OK) CxxErrorHandler();
        vbat += HAL_ADC_GetValue(&hadc1);
    }

    vbat /= 8;

    if (HAL_ADC_Stop(&hadc1) != HAL_OK) CxxErrorHandler();
    if (HAL_TIM_Base_Stop(&htim6) != HAL_OK)
        CxxErrorHandler();

    htim6.Init.Period = 1817;
    if (HAL_TIM_Base_Init(&htim6) != HAL_OK) CxxErrorHandler();

    HAL_Delay(1);

    gpio::BAT_DIVIDER::on();
    if (!usb_ce) gpio::USB_CE::off();   // Restore battery charging state.

    INFO("Vref = %lu", vrefint);
    INFO("Vbat = %lu (raw)", vbat);

    // Order of operations is important to avoid underflow.
    vbat *= 6600;
    vbat /= (vref + 1);

    INFO("Vbat = %lumV", vbat);

    uint8_t data[3];
    data[0] = kiss::hardware::GET_BATTERY_LEVEL;
    data[1] = (vbat >> 8) & 0xFF;
    data[2] = (vbat & 0xFF);

    ioport->write(data, 3, 6, 10);
    DEBUG("exit pollBatteryLevel");
}

#if 0
void stop() {
    osDelay(100);
#if 0
    auto restore = SysTick->CTRL;

    kiss::settings().input_offset += 6;
    setAudioInputLevels();
    kiss::settings().input_offset -= 6;
    DEBUG("Stop");
    // __disable_irq();
    vTaskSuspendAll();
    SysTick->CTRL = 0;
    HAL_COMP_Init(&hcomp1);
    HAL_COMP_Start_IT(&hcomp1);
    while (adcState == STOPPED) {
        // PWR_MAINREGULATOR_ON / PWR_LOWPOWERREGULATOR_ON
        HAL_PWR_EnterSTOPMode(PWR_LOWPOWERREGULATOR_ON, PWR_STOPENTRY_WFI);
    }
    SystemClock_Config();
    SysTick->CTRL = restore;
    // __enable_irq();
    HAL_COMP_Stop_IT(&hcomp1);
    HAL_COMP_DeInit(&hcomp1);
    xTaskResumeAll();
    setAudioInputLevels();
    // adcState = DEMODULATOR;
    DEBUG("Wake");
#endif
}
#endif

}}} // mobilinkd::tnc::audio
