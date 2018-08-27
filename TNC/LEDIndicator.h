// Copyright 2017 Rob Riggs <rob@mobilinkd.com>
// All rights reserved.

#ifndef MOBILINKD__TNC__LED_INDICATOR_H_
#define MOBILINKD__TNC__LED_INDICATOR_H_

#ifdef __cplusplus
extern "C" {
#endif

void indicate_turning_on(void);
void indicate_on(void);
void indicate_initializing_ble(void);
void indicate_turning_off(void);
void indicate_waiting_to_connect(void);
void indicate_connected_via_usb(void);
void indicate_connected_via_ble(void);
void tx_on(void);
void tx_off(void);
void rx_on(void);
void rx_off(void);

void HTIM1_PeriodElapsedCallback(void);

#ifdef __cplusplus
}

#endif //__cplusplus

#endif // MOBILINKD__TNC__LED_INDICATOR_H_
