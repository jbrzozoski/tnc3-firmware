// Copyright 2020 Rob Riggs <rob@mobilinkd.com>
// All rights reserved.

#include "Demodulator.hpp"

namespace mobilinkd { namespace tnc {

/**
 * Start the ADC DMA transfer.  The block size is equal to the number of
 * 32-bit elements in the buffer.  This is also equal to the number of
 * 16-bit elements in each DMA "half complete" transfer.  Each DMA transfer
 * results in block_size * 2 bytes being transferred.
 *
 * We must bear in mind that the DMA buffer size is expressed in DWORDs,
 * the DMA transfer size is expressed in WORDs, and the memory copy operation
 * in BYTEs.
 *
 * @param period
 * @param block_size
 */
void IDemodulator::startADC(uint32_t period, uint32_t block_size)
{
    audio::set_adc_block_size(block_size);

    htim6.Init.Period = period;
    htim6.Init.AutoReloadPreload = TIM_AUTORELOAD_PRELOAD_DISABLE;
    if (HAL_TIM_Base_Init(&htim6) != HAL_OK)
    {
        CxxErrorHandler();
    }

    if (HAL_TIM_Base_Start(&htim6) != HAL_OK)
    {
        CxxErrorHandler();
    }

    if (HAL_ADC_Start_DMA(&hadc1, audio::adc_buffer,
        block_size * 2) != HAL_OK)
    {
        CxxErrorHandler();
    }
}

}} // mobilinkd::tnc