// Copyright 2015-2019 Mobilinkd LLC <rob@mobilinkd.com>
// All rights reserved.


#ifndef MOBILINKD__HDLC_FRAME_HPP_
#define MOBILINKD__HDLC_FRAME_HPP_

#ifndef EXCLUDE_CRC
#include "stm32l4xx_hal.h"
extern CRC_HandleTypeDef hcrc;
#endif

#include "cmsis_os.h"

#include <Log.h>
#include "SegmentedBuffer.hpp"

#include <boost/intrusive/list.hpp>

#include <iterator>
#include <algorithm>

namespace mobilinkd { namespace tnc { namespace hdlc {


using boost::intrusive::list_base_hook;
using boost::intrusive::list;
using boost::intrusive::constant_time_size;

template <typename POOL, POOL* allocator>
class Frame : public list_base_hook<>
{
public:
    typedef POOL pool_type;
    typedef buffer::SegmentedBuffer<POOL, allocator> data_type;
    typedef typename data_type::iterator iterator;

    enum Type {
        DATA = 0, TX_DELAY, P_PERSIST, SLOT_TIME, TX_TAIL, DUPLEX, HARDWARE,
        TEXT, LOG};

    enum Source {
      RF_DATA = 0x00, SERIAL_DATA = 0x10, DIGI_DATA = 0x20,
      FRAME_RETURN = 0xF0};

private:
    data_type data_;
    int crc_{-1};
    int fcs_{-2};
    bool complete_{false};
    uint8_t frame_type_{Type::DATA};

#ifndef EXCLUDE_CRC
    uint16_t compute_crc(iterator first) {

        uint16_t bytes = data_.size();
        uint16_t block_size = std::min(bytes, POOL::chunk_type::size());

        uint32_t checksum = HAL_CRC_Calculate(&hcrc, (uint32_t*) &(*first), block_size);

        bytes -= block_size;
        std::advance(first, block_size);

        while (bytes) {
            block_size = std::min(bytes, POOL::chunk_type::size());
            checksum = HAL_CRC_Accumulate(&hcrc, (uint32_t*) &(*first), block_size);
            bytes -= block_size;
            std::advance(first, block_size);
        }

        checksum ^= 0xFFFF;  // Compliment
        checksum <<= 16;     // Shift
        asm volatile("rbit %0, %0" : "+r" (checksum)); // Reverse
        uint16_t result = checksum & 0xFFFF;
        DEBUG("CRC = %hx", result);
        return result;
    }
#else
    uint16_t compute_crc(iterator first, iterator last) {return 0;}
#endif

public:
    Frame()
    : list_base_hook<>(), data_(), crc_(-1), fcs_(-2), complete_(false)
    {}

    uint8_t type() const {return frame_type_ & 0x0F;}
    void type(uint8_t t) {frame_type_ |= ((frame_type_ & 0xF0) | t);}

    uint8_t source() const {return frame_type_ & 0xF0;}
    void source(uint8_t s) {frame_type_ |= ((frame_type_ & 0x0F) | s);}

    void clear() {
        data_.clear();
        crc_ = -1;
        fcs_ = -2;
        complete_ = false;
        frame_type_ = 0;    // RF_DATA.
    }

    void assign(data_type& data) {
        data_.splice(data_.end(), data);
    }

    uint16_t size() const {return data_.size();}

    uint16_t crc() const {return crc_;}
    uint16_t fcs() const {return fcs_;}

    bool complete() const {return complete_;}

    bool ok() const {return crc_ == 0x0f47; /*0xf0b8;*/}

    typename data_type::iterator begin() { return data_.begin(); }
    typename data_type::iterator end() { return data_.end(); }

    bool push_back(uint8_t value)
    {
        return data_.push_back(value);
    }

    void add_fcs() {           // TX frames have the checksums added.
        fcs_ = compute_crc(data_.begin());
        data_.push_back(uint8_t(fcs_ & 0xFF));
        data_.push_back(uint8_t((fcs_ >> 8) & 0xFF));
        crc_ = 0x0f47;
        complete_ = true;
    }

    void parse_fcs() {              // RX frames have the checksums parsed.
        auto it = data_.begin();
        std::advance(it, data_.size() - 2);
        fcs_ = (*it);
        ++it;
        fcs_ |= (*it) << 8;
        DEBUG("FCS = %hx", fcs_);
        crc_ = compute_crc(data_.begin());
        complete_ = true;
    }

    void frame_and_stuff_data() {
        // I went back and forth several times on how to do this "cleanly"...
        // - I initially wanted to just allocate a whole new IOFrame up in the caller, but then realized I couldn't duplicate the meta-data like crc, fcs, frame_type, complete...
        // - Then I wanted to replace the contents of the current data_, but it's not a pointer and I can't just hijack it...
        // - So, for the time being, I'm just double-copying out to a temp data_type and then back in.  It sucks, but I wanted to get something done and prove it worked before working to make it efficient...
        data_type temp_buffer;

        // TODO !!! Add FLAG bytes (0x7E) at beginning of frame
        temp_buffer.push_back(0x7E);

        // TODO !!! Walk full size of data, and do bit stuffing
        uint8_t output_byte = 0x00;
        int output_bits = 0;
        auto bytes_to_read = data_.size();
        auto data_in = data_.begin();
        int consecutive_ones = 0;
        uint8_t byte;
        while (bytes_to_read) {
            byte = *data_in;
            std::advance(data_in, 1);
            bytes_to_read--;
            for (int i = 0; i != 8; i++) {
                uint8_t bit = byte & 1;

                if (bit) {
                    output_byte |= (0x01 << output_bits);
                    ++consecutive_ones;
                } else {
                    consecutive_ones = 0;
                }
                ++output_bits;
                if (output_bits >= 8) {
                    temp_buffer.push_back(output_byte);
                    output_byte = 0x00;
                    output_bits = 0;
                }
                if (consecutive_ones >= 5) {
                    ++output_bits;
                    if (output_bits >= 8) {
                        temp_buffer.push_back(output_byte);
                        output_byte = 0x00;
                        output_bits = 0;
                    }
                    consecutive_ones = 0;
                }
                byte >>= 1;
            }
        }

        // Make sure we have one FLAG byte in matching bit-stuffing alignment at the end of the frame
        byte = 0x7E;
        for (int i = 0; i != 8; i++) {
            uint8_t bit = byte & 1;

            if (bit) {
                output_byte |= (0x01 << output_bits);
            }
            ++output_bits;
            if (output_bits >= 8) {
                temp_buffer.push_back(output_byte);
                output_byte = 0x00;
                output_bits = 0;
            }
            byte >>= 1;
        }

        // Add partial FLAG bits to pad to a full byte...
        byte = 0x7E;
        while (output_bits != 0)
        {
            uint8_t bit = byte & 1;

            if (bit) {
                output_byte |= (0x01 << output_bits);
            }
            ++output_bits;
            if (output_bits >= 8) {
                temp_buffer.push_back(output_byte);
                output_byte = 0x00;
                output_bits = 0;
            }
            byte >>= 1;
        }

        // TODO - Add extra FLAG and 0 padding to the over-the-air and RS sizes if we are going to be doing FX.25?

        // TODO - Now copy this data back into the IOFrame buffer...
        data_.clear();
        for (auto c : temp_buffer) data_.push_back(c);
    }

};

template <typename Frame, size_t SIZE = 16>
class FramePool
{
public:
    typedef Frame frame_type;
    typedef list<frame_type, constant_time_size<false>> FrameList;

private:
    static const uint16_t FRAME_COUNT = SIZE;
    frame_type frames_[FRAME_COUNT];
    FrameList free_list_;

public:
    FramePool()
    : frames_(), free_list_()
    {
        for (auto& frame : frames_) {
            free_list_.push_back(frame);
        }
    }

    uint16_t size() const {return free_list_.size();}

    static constexpr uint16_t capacity() { return SIZE; }

    frame_type* acquire() {
        frame_type* result = nullptr;
        auto x = taskENTER_CRITICAL_FROM_ISR();
        if (!free_list_.empty()) {
            result = &free_list_.front();
            free_list_.pop_front();
        }
        taskEXIT_CRITICAL_FROM_ISR(x);
        DEBUG("Acquired frame %p (size after = %d)", result, free_list_.size());
        return result;
    }

    void release(frame_type* frame) {
        DEBUG("Released frame %p (size before = %d)", frame, free_list_.size());
        frame->clear();
        auto x = taskENTER_CRITICAL_FROM_ISR();
        free_list_.push_back(*frame);
        taskEXIT_CRITICAL_FROM_ISR(x);
    }
};

typedef buffer::Pool<48> FrameSegmentPool;  // 12K buffer of frames;

extern FrameSegmentPool frameSegmentPool;

typedef Frame<FrameSegmentPool, &frameSegmentPool> IoFrame;
typedef FramePool<IoFrame, 48> IoFramePool;

IoFramePool& ioFramePool(void);

/**
 * This is needed to work around a major defect in the GCC 5.2.0 compiler
 * that causes functions using ioFramePool().release(frame) to use an
 * extremely large amount of stack space (4-6K vs. 24 bytes).
 *
 * This function merely acts as a compiler firewall to the stack usage.
 * It's own stack usage is minimal even though it is making the same call.
 *
 * @param frame
 */
void release(IoFrame* frame);

IoFrame* acquire(void);
IoFrame* acquire_wait(void);

}}} // mobilinkd::tnc::hdlc

#endif // MOBILINKD__HDLC_FRAME_HPP_
