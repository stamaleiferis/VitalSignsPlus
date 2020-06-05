#ifndef COEFFS_H
#define COEFFS_H

float NOTCH_60_A[] = {0.60380455, 0.95395256};
float NOTCH_60_B[] =  {0.97697628, 0.60380455, 0.97697628};
float NOTCH_80_A[] =  {1.56873452, 0.93906251};
float NOTCH_80_B[] = {0.96953125, 1.56873452, 0.96953125};



float TAPS_KAISER10[] = {-0.0010072901142833142, 0.003601265264628422, 0.008345248716088195,
0.006317226962764452,
-0.00013749849962718056,
-0.0024050901939474626,
0.002513211553054252,
0.007595518430484113,
0.005313900641389956,
-0.001803527147671388,
-0.004360626936739569,
0.0009484692716605278,
0.006459928115682073,
0.003902895285244263,
-0.003981741654643091,
-0.006859511298095418,
-0.001048308142361605,
0.005017065400334756,
0.0021602687812766917,
-0.006627449114542474,
-0.00986544347867477,
-0.0034002273953086556,
0.0033903078637572627,
0.00020200396858800502,
-0.009675752689803767,
-0.013330721975049998,
-0.006002108774380662,
0.0017524088289586175,
-0.0018145435295529432,
-0.01305408058608369,
-0.017215325018074434,
-0.008726651855489691,
0.00033651812730625674,
-0.003685042417069135,
-0.01670402815418874,
-0.021521917363994357,
-0.01143239866472341,
-0.000534585818319052,
-0.005139410302440773,
-0.02062247658072226,
-0.02636714498289277,
-0.013972931701328668,
-0.0003680586958750374,
-0.005784054489420107,
-0.024953308063387903,
-0.032154651994723144,
-0.016206639160195088,
0.0017423027155777163,
-0.004914197851898864,
-0.03024910405617826,
-0.04011099708559525,
-0.01800632632314211,
0.008019877119400068,
-0.000761573117313688,
-0.03852784766185687,
-0.05467807136077638,
-0.019267955240082305,
0.02712532178310358,
0.014413495985148232,
-0.061311885326081075,
-0.10633366908484555,
-0.019917854704871212,
0.18185005663659862,
0.3541990032212447,
0.3541990032212447,
0.18185005663659862,
-0.019917854704871212,
-0.10633366908484555,
-0.061311885326081075,
0.014413495985148232,
0.02712532178310358,
-0.019267955240082305,
-0.05467807136077638,
-0.03852784766185687,
-0.000761573117313688,
0.008019877119400068,
-0.01800632632314211,
-0.04011099708559525,
-0.03024910405617826,
-0.004914197851898864,
0.0017423027155777163,
-0.016206639160195088,
-0.032154651994723144,
-0.024953308063387903,
-0.005784054489420107,
-0.0003680586958750374,
-0.013972931701328668,
-0.02636714498289277,
-0.02062247658072226,
-0.005139410302440773,
-0.000534585818319052,
-0.01143239866472341,
-0.021521917363994357,
-0.01670402815418874,
-0.003685042417069135,
0.00033651812730625674,
-0.008726651855489691,
-0.017215325018074434,
-0.01305408058608369,
-0.0018145435295529432,
0.0017524088289586175,
-0.006002108774380662,
-0.013330721975049998,
-0.009675752689803767,
0.00020200396858800502,
0.0033903078637572627,
-0.0034002273953086556,
-0.00986544347867477,
-0.006627449114542474,
0.0021602687812766917,
0.005017065400334756,
-0.001048308142361605,
-0.006859511298095418,
-0.003981741654643091,
0.003902895285244263,
0.006459928115682073,
0.0009484692716605278,
-0.004360626936739569,
-0.001803527147671388,
0.005313900641389956,
0.007595518430484113,
0.002513211553054252,
-0.0024050901939474626,
-0.00013749849962718056,
0.006317226962764452,
0.008345248716088195,
0.003601265264628422,
-0.0010072901142833142};
float TAPS_KAISER16[] =  {-0.0008586890456536407,
0.003085935109673783,
0.007187533823102203,
0.005468081239621113,
-0.0001196006447588263,
-0.002102100755336547,
0.0022069759824999563,
0.006700889511430655,
0.004709288390965654,
-0.0016054322075336148,
-0.0038985805711122664,
0.0008515894071998858,
0.00582432267224028,
0.0035332712860013326,
-0.0036190728054081916,
-0.006259124774552797,
-0.0009602135601027526,
0.004612639772965845,
0.001993384742990657,
-0.006137283206916916,
-0.009167604040960155,
-0.0031704426693012554,
0.003171664214535173,
0.00018958661813879194,
-0.00910953080965035,
-0.012589031350267241,
-0.005685045984253635,
0.001664643593682293,
-0.0017285156538818871,
-0.012469149656484997,
-0.016487509656477554,
-0.008379181077082272,
0.0003239227177954395,
-0.003555654128697858,
-0.016155019245268246,
-0.020861319068800028,
-0.011105499798237679,
-0.000520383638019536,
-0.005012917462253873,
-0.020153703125377644,
-0.02581545148284545,
-0.013704797692334946,
-0.00036160683626784473,
-0.00569183687902235,
-0.024593189703435214,
-0.031736821035966076,
-0.01601812710652785,
0.0017242804421477536,
-0.004869326060970294,
-0.030007315090967823,
-0.03983299962509072,
-0.017899305534926004,
0.007979522722731667,
-0.0007583778042033109,
-0.03839549513921852,
-0.054527674197619884,
-0.019226686957435995,
0.027081678161390055,
0.014396889944140002,
-0.06126459366002082,
-0.10628404910081365,
-0.01991311230226165,
0.18183446867920314,
0.35419562965456486,
0.35419562965456486,
0.18183446867920314,
-0.01991311230226165,
-0.10628404910081365,
-0.06126459366002082,
0.014396889944140002,
0.027081678161390055,
-0.019226686957435995,
-0.054527674197619884,
-0.03839549513921852,
-0.0007583778042033109,
0.007979522722731667,
-0.017899305534926004,
-0.03983299962509072,
-0.030007315090967823,
-0.004869326060970294,
0.0017242804421477536,
-0.01601812710652785,
-0.031736821035966076,
-0.024593189703435214,
-0.00569183687902235,
-0.00036160683626784473,
-0.013704797692334946,
-0.02581545148284545,
-0.020153703125377644,
-0.005012917462253873,
-0.000520383638019536,
-0.011105499798237679,
-0.020861319068800028,
-0.016155019245268246,
-0.003555654128697858,
0.0003239227177954395,
-0.008379181077082272,
-0.016487509656477554,
-0.012469149656484997,
-0.0017285156538818871,
0.001664643593682293,
-0.005685045984253635,
-0.012589031350267241,
-0.00910953080965035,
0.00018958661813879194,
0.003171664214535173,
-0.0031704426693012554,
-0.009167604040960155,
-0.006137283206916916,
0.001993384742990657,
0.004612639772965845,
-0.0009602135601027526,
-0.006259124774552797,
-0.0036190728054081916,
0.0035332712860013326,
0.00582432267224028,
0.0008515894071998858,
-0.0038985805711122664,
-0.0016054322075336148,
0.004709288390965654,
0.006700889511430655,
0.0022069759824999563,
-0.002102100755336547,
-0.0001196006447588263,
0.005468081239621113,
0.007187533823102203,
0.003085935109673783,
-0.0008586890456536407};

#endif
