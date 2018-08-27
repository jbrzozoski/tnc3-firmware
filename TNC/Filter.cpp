#include "Filter.h"

#include <cstdlib>

const float filter_taps[13][FILTER_TAP_NUM] = {
    {
        // 0dB Hamming (1693.3)
        0.0092140576627,
        0.00768806806252,
        0.00498340794238,
        -0.00171692105003,
        -0.0142537522489,
        -0.0319314833665,
        -0.0509364246735,
        -0.0650918044189,
        -0.0677825846114,
        -0.0544018818,
        -0.0244193735639,
        0.0177467007765,
        0.0634791554618,
        0.102191070678,
        0.124356477862,
        0.124356477862,
        0.102191070678,
        0.0634791554618,
        0.0177467007765,
        -0.0244193735639,
        -0.0544018818,
        -0.0677825846114,
        -0.0650918044189,
        -0.0509364246735,
        -0.0319314833665,
        -0.0142537522489,
        -0.00171692105003,
        0.00498340794238,
        0.00768806806252,
        0.0092140576627,
    },
    {
        // 1dB Hamming (1749.2)
        0.00984360063296,
        0.00907660156831,
        0.00733802154094,
        0.00165885599255,
        -0.01027803496,
        -0.028288428514,
        -0.0487844743735,
        -0.0653441331879,
        -0.0706814843758,
        -0.0593532384181,
        -0.0301832127721,
        0.0125781304703,
        0.0599187881831,
        0.100456125433,
        0.123807723158,
        0.123807723158,
        0.100456125433,
        0.0599187881831,
        0.0125781304703,
        -0.0301832127721,
        -0.0593532384181,
        -0.0706814843758,
        -0.0653441331879,
        -0.0487844743735,
        -0.028288428514,
        -0.01027803496,
        0.00165885599255,
        0.00733802154094,
        0.00907660156831,
        0.00984360063296,
    },
    {
        // 2dB Hamming (1805.6)
        0.0101071938454,
        0.0101751907488,
        0.00949571140771,
        0.00500536790767,
        -0.00608985117047,
        -0.0241824665512,
        -0.0460028160886,
        -0.0649487308381,
        -0.0730641303456,
        -0.0640057046058,
        -0.0358585166096,
        0.00736846731409,
        0.0562886910258,
        0.098688560769,
        0.12326772148,
        0.12326772148,
        0.098688560769,
        0.0562886910258,
        0.00736846731409,
        -0.0358585166096,
        -0.0640057046058,
        -0.0730641303456,
        -0.0649487308381,
        -0.0460028160886,
        -0.0241824665512,
        -0.00608985117047,
        0.00500536790767,
        0.00949571140771,
        0.0101751907488,
        0.0101071938454,
    },
    {
        // 3dB Hamming (1861.5)
        0.00999732596792,
        0.0109364331138,
        0.0113663405297,
        0.00819091080479,
        -0.00184492788726,
        -0.0197616766146,
        -0.0426984518018,
        -0.0639498932619,
        -0.0749094986595,
        -0.0682876192213,
        -0.0413467222516,
        0.00222026614421,
        0.0526805390514,
        0.0969649735461,
        0.1228029701,
        0.1228029701,
        0.0969649735461,
        0.0526805390514,
        0.00222026614421,
        -0.0413467222516,
        -0.0682876192213,
        -0.0749094986595,
        -0.0639498932619,
        -0.0426984518018,
        -0.0197616766146,
        -0.00184492788726,
        0.00819091080479,
        0.0113663405297,
        0.0109364331138,
        0.00999732596792,
    },
    {
        // 4dB Hamming (1917.0)
        0.00953170942079,
        0.0113523463692,
        0.012917947221,
        0.0111611200433,
        0.00238829862135,
        -0.0150956296368,
        -0.0389313364275,
        -0.0623942365343,
        -0.07625339099,
        -0.0722276074899,
        -0.0466690276902,
        -0.00287524097904,
        0.0491038547992,
        0.095314261046,
        0.122454987694,
        0.122454987694,
        0.095314261046,
        0.0491038547992,
        -0.00287524097904,
        -0.0466690276902,
        -0.0722276074899,
        -0.07625339099,
        -0.0623942365343,
        -0.0389313364275,
        -0.0150956296368,
        0.00238829862135,
        0.0111611200433,
        0.012917947221,
        0.0113523463692,
        0.00953170942079,
    },
    {
        // 5dB Hamming (1972.5)
        0.00872938144139,
        0.011420928511,
        0.0141310255503,
        0.0138829549523,
        0.00657217070647,
        -0.0102167866884,
        -0.0347251133165,
        -0.0603020097307,
        -0.0771226785494,
        -0.0758651302103,
        -0.0518731312861,
        -0.00796039989017,
        0.0455363026357,
        0.0937399215189,
        0.122245216908,
        0.122245216908,
        0.0937399215189,
        0.0455363026357,
        -0.00796039989017,
        -0.0518731312861,
        -0.0758651302103,
        -0.0771226785494,
        -0.0603020097307,
        -0.0347251133165,
        -0.0102167866884,
        0.00657217070647,
        0.0138829549523,
        0.0141310255503,
        0.011420928511,
        0.00872938144139,
    },
    {
        // 6dB Hamming (2027.0)
        0.00763984171566,
        0.0111488308293,
        0.0149655370284,
        0.0162631703876,
        0.0105693777806,
        -0.00527673382432,
        -0.0302054840484,
        -0.0577401237325,
        -0.077513860058,
        -0.079139201575,
        -0.0568701896276,
        -0.0129510233888,
        0.0420349410452,
        0.0922665002,
        0.122176947883,
        0.122176947883,
        0.0922665002,
        0.0420349410452,
        -0.0129510233888,
        -0.0568701896276,
        -0.079139201575,
        -0.077513860058,
        -0.0577401237325,
        -0.0302054840484,
        -0.00527673382432,
        0.0105693777806,
        0.0162631703876,
        0.0149655370284,
        0.0111488308293,
        0.00763984171566,
    },
    {
        // 7dB Cosine (1986.0)
        0.00460770532874,
        0.0162873765157,
        0.0251984580718,
        0.0248176615108,
        0.0117845748105,
        -0.0125423969002,
        -0.0421426113528,
        -0.0680353327511,
        -0.0811090126152,
        -0.0751714219734,
        -0.0492188545506,
        -0.00815365116672,
        0.038335087686,
        0.0784276419333,
        0.101570441881,
        0.101570441881,
        0.0784276419333,
        0.038335087686,
        -0.00815365116672,
        -0.0492188545506,
        -0.0751714219734,
        -0.0811090126152,
        -0.0680353327511,
        -0.0421426113528,
        -0.0125423969002,
        0.0117845748105,
        0.0248176615108,
        0.0251984580718,
        0.0162873765157,
        0.00460770532874,
    },
    {
        // 8dB Cosine (2021.0, 2021.1)
        0.00420597097315,
        0.0159601020805,
        0.0260190570929,
        0.0273151975444,
        0.0157146317097,
        -0.00808663238039,
        -0.0383509470879,
        -0.0659043862717,
        -0.0810670338405,
        -0.0769492259504,
        -0.0520051674945,
        -0.0109566710433,
        0.0362780863637,
        0.0773640302855,
        0.101184010968,
        0.101184010968,
        0.0773640302855,
        0.0362780863637,
        -0.0109566710433,
        -0.0520051674945,
        -0.0769492259504,
        -0.0810670338405,
        -0.0659043862717,
        -0.0383509470879,
        -0.00808663238039,
        0.0157146317097,
        0.0273151975444,
        0.0260190570929,
        0.0159601020805,
        0.00420597097315,

    },
    {
        // 9dB Hamming (2186.5)
        0.00308518474955,
        0.00850291315071,
        0.0151934174423,
        0.0211164164746,
        0.0210280537958,
        0.00946135065651,
        -0.0151061442669,
        -0.0474011160773,
        -0.0758516294354,
        -0.0867915071262,
        -0.0707058729109,
        -0.0275514181864,
        0.0316531744803,
        0.0881474710462,
        0.122515729369,
        0.122515729369,
        0.0881474710462,
        0.0316531744803,
        -0.0275514181864,
        -0.0707058729109,
        -0.0867915071262,
        -0.0758516294354,
        -0.0474011160773,
        -0.0151061442669,
        0.00946135065651,
        0.0210280537958,
        0.0211164164746,
        0.0151934174423,
        0.00850291315071,
        0.00308518474955,
    },
    {
        // 10dB Cosine (2087.0, 2087.1)
        0.003297704105,
        0.0148243949961,
        0.0268172892254,
        0.0313287864329,
        0.0227440831014,
        0.00036329489673,
        -0.0307896841323,
        -0.061300901507,
        -0.0804397189579,
        -0.0799281940065,
        -0.0570830747612,
        -0.0161917885552,
        0.0324263666501,
        0.0754301963068,
        0.100574174183,
        0.100574174183,
        0.0754301963068,
        0.0324263666501,
        -0.0161917885552,
        -0.0570830747612,
        -0.0799281940065,
        -0.0804397189579,
        -0.061300901507,
        -0.0307896841323,
        0.00036329489673,
        0.0227440831014,
        0.0313287864329,
        0.0268172892254,
        0.0148243949961,
        0.003297704105,
    },
    {
        // 11dB Cosine (2119.0, 2119.1)
        0.00279757679371,
        0.0140443163212,
        0.0268514195956,
        0.0329271588621,
        0.0259396342176,
        0.00444730813715,
        -0.0269594711541,
        -0.0588146601004,
        -0.0798925711826,
        -0.0812085944654,
        -0.0594715656311,
        -0.01871136648,
        0.0305735278948,
        0.0745360339334,
        0.100349256713,
        0.100349256713,
        0.0745360339334,
        0.0305735278948,
        -0.01871136648,
        -0.0594715656311,
        -0.0812085944654,
        -0.0798925711826,
        -0.0588146601004,
        -0.0269594711541,
        0.00444730813715,
        0.0259396342176,
        0.0329271588621,
        0.0268514195956,
        0.0140443163212,
        0.00279757679371,
    },
    {
        // 12dB Cosine (2150,2151)
        0.0022752196356,
        0.0131417834537,
        0.0266616639804,
        0.0342658094422,
        0.028925742996,
        0.00843279814311,
        -0.0231061648363,
        -0.0562188716332,
        -0.0792056020955,
        -0.0823667306435,
        -0.0617763153975,
        -0.0211791883525,
        0.0287608072456,
        0.0736871278906,
        0.100177853829,
        0.100177853829,
        0.0736871278906,
        0.0287608072456,
        -0.0211791883525,
        -0.0617763153975,
        -0.0823667306435,
        -0.0792056020955,
        -0.0562188716332,
        -0.0231061648363,
        0.00843279814311,
        0.028925742996,
        0.0342658094422,
        0.0266616639804,
        0.0131417834537,
        0.0022752196356,
    }
};

void Filter_init(Filter* f, const float* taps) {
    for (size_t i = 0; i != FILTER_TAP_NUM; ++i) {
        f->history[i] = 0;
        f->taps[i] = taps[i];   // initialize ccmram.
    }
    f->last_index = 0;
}

void Filter_put(Filter* f, float input) {
    f->history[f->last_index++] = input;
    if (f->last_index == FILTER_TAP_NUM)
        f->last_index = 0;
}

float Filter_get(Filter* f) {
    float acc = 0;
    int index = f->last_index, i;
    for (i = 0; i < FILTER_TAP_NUM; ++i) {
        index = index != 0 ? index - 1 : FILTER_TAP_NUM - 1;
        acc += f->history[index] * f->taps[i];
    };
    return acc;
}
