// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#include "h2d_common.h"
#include "shapeset_common.h"
#include "shapeset_h1_all.h"

//eigen functions for Laplace

// p = 2
#define eigen_laplace_p2_2(x)   ((1.0000000000000000) * l2(x) )


// p = 3
#define eigen_laplace_p3_2(x)   ((1.0000000000000000) * l3(x) )
#define eigen_laplace_p3_3(x)   ((1.0000000000000000) * l2(x) )


// p = 4
#define eigen_laplace_p4_2(x)   ((0.1200767132565378) * l2(x) + (0.9927646160764934) * l4(x) )
#define eigen_laplace_p4_3(x)   ((-0.9999999999999999) * l3(x) )
#define eigen_laplace_p4_4(x)   ((0.9927646160764932) * l2(x) + (-0.1200767132565379) * l4(x) )


// p = 5
#define eigen_laplace_p5_2(x)   ((-0.2721806201771906) * l3(x) + (-0.9622461795195447) * l5(x) )
#define eigen_laplace_p5_3(x)   ((-0.1200767132565378) * l2(x) + (-0.9927646160764934) * l4(x) )
#define eigen_laplace_p5_4(x)   ((-0.9622461795195447) * l3(x) + (0.2721806201771906) * l5(x) )
#define eigen_laplace_p5_5(x)   ((-0.9927646160764934) * l2(x) + (0.1200767132565378) * l4(x) )


// p = 6
#define eigen_laplace_p6_2(x)   ((0.0460292214008892) * l2(x) + (0.4098443624284846) * l4(x) + (0.9109933640608021) * l6(x) )
#define eigen_laplace_p6_3(x)   ((0.2721806201771907) * l3(x) + (0.9622461795195447) * l5(x) )
#define eigen_laplace_p6_4(x)   ((0.1111170386828487) * l2(x) + (0.9041998053176488) * l4(x) + (-0.4124023712078773) * l6(x) )
#define eigen_laplace_p6_5(x)   ((0.9622461795195448) * l3(x) + (-0.2721806201771906) * l5(x) )
#define eigen_laplace_p6_6(x)   ((0.9927408093211348) * l2(x) + (-0.1202094449247417) * l4(x) + (0.0039210788443064) * l6(x) )


// p = 7
#define eigen_laplace_p7_2(x)   ((0.1244273041355310) * l3(x) + (0.5154918741058946) * l5(x) + (0.8478124637656356) * l7(x) )
#define eigen_laplace_p7_3(x)   ((-0.0460292214008892) * l2(x) + (-0.4098443624284848) * l4(x) + (-0.9109933640608022) * l6(x) )
#define eigen_laplace_p7_4(x)   ((0.2461133039997421) * l3(x) + (0.8117163656110983) * l5(x) + (-0.5296647839845884) * l7(x) )
#define eigen_laplace_p7_5(x)   ((0.1111170386828487) * l2(x) + (0.9041998053176490) * l4(x) + (-0.4124023712078774) * l6(x) )
#define eigen_laplace_p7_6(x)   ((0.9612211439517415) * l3(x) + (-0.2745626877962528) * l5(x) + (0.0258697292255408) * l7(x) )
#define eigen_laplace_p7_7(x)   ((-0.9927408093211348) * l2(x) + (0.1202094449247416) * l4(x) + (-0.0039210788443064) * l6(x) )


// p = 8
#define eigen_laplace_p8_2(x)   ((0.0233791285597233) * l2(x) + (0.2118352349060833) * l4(x) + (0.5868165705262892) * l6(x) + (0.7811693556174851) * l8(x) )
#define eigen_laplace_p8_3(x)   ((0.1244273041355311) * l3(x) + (0.5154918741058940) * l5(x) + (0.8478124637656346) * l7(x) )
#define eigen_laplace_p8_4(x)   ((0.0418264255244555) * l2(x) + (0.3682779654706578) * l4(x) + (0.6912121989219027) * l6(x) + (-0.6203608517130881) * l8(x) )
#define eigen_laplace_p8_5(x)   ((0.2461133039997420) * l3(x) + (0.8117163656110977) * l5(x) + (-0.5296647839845874) * l7(x) )
#define eigen_laplace_p8_6(x)   ((-0.1103162273830400) * l2(x) + (-0.8972440333626744) * l4(x) + (0.4217305141658119) * l6(x) + (-0.0701915094045806) * l8(x) )
#define eigen_laplace_p8_7(x)   ((-0.9612211439517412) * l3(x) + (0.2745626877962528) * l5(x) + (-0.0258697292255406) * l7(x) )
#define eigen_laplace_p8_8(x)   ((-0.9927408002351916) * l2(x) + (0.1202094754840794) * l4(x) + (-0.0039219942315734) * l6(x) + (0.0000592746276296) * l8(x) )


// p = 9
#define eigen_laplace_p9_2(x)   ((-0.0679631762828228) * l3(x) + (-0.2941064986189726) * l5(x) + (-0.6288879599997316) * l7(x) + (-0.7165070187423248) * l9(x) )
#define eigen_laplace_p9_3(x)   ((-0.0233791285597233) * l2(x) + (-0.2118352349060834) * l4(x) + (-0.5868165705262897) * l6(x) + (-0.7811693556174855) * l8(x) )
#define eigen_laplace_p9_4(x)   ((-0.1166010439099156) * l3(x) + (-0.4632814174552214) * l5(x) + (-0.5506359890315790) * l7(x) + (0.6845250414579539) * l9(x) )
#define eigen_laplace_p9_5(x)   ((0.0418264255244555) * l2(x) + (0.3682779654706578) * l4(x) + (0.6912121989219027) * l6(x) + (-0.6203608517130882) * l8(x) )
#define eigen_laplace_p9_6(x)   ((0.2405135991654621) * l3(x) + (0.7896100977856557) * l5(x) + (-0.5482950511284391) * l7(x) + (0.1343191683992437) * l9(x) )
#define eigen_laplace_p9_7(x)   ((0.1103162273830400) * l2(x) + (0.8972440333626747) * l4(x) + (-0.4217305141658118) * l6(x) + (0.0701915094045807) * l8(x) )
#define eigen_laplace_p9_8(x)   ((0.9612171512437465) * l3(x) + (-0.2745680046649571) * l5(x) + (0.0259321857669741) * l7(x) + (-0.0012331709159256) * l9(x) )
#define eigen_laplace_p9_9(x)   ((-0.9927408002351915) * l2(x) + (0.1202094754840794) * l4(x) + (-0.0039219942315734) * l6(x) + (0.0000592746276296) * l8(x) )


// p = 10
#define eigen_laplace_p10_2(x)   ((-0.0135012819178048) * l2(x) + (-0.1230733557674822) * l4(x) + (-0.3644642014036448) * l6(x) + (-0.6487214117012139) * l8(x) + (-0.6565036484150675) * l10(x) )
#define eigen_laplace_p10_3(x)   ((-0.0679631762828231) * l3(x) + (-0.2941064986189725) * l5(x) + (-0.6288879599997312) * l7(x) + (-0.7165070187423248) * l9(x) )
#define eigen_laplace_p10_4(x)   ((-0.0229617482557899) * l2(x) + (-0.2064896698810138) * l4(x) + (-0.5229445560893998) * l6(x) + (-0.3991433171533719) * l8(x) + (0.7239118578789349) * l10(x) )
#define eigen_laplace_p10_5(x)   ((-0.1166010439099159) * l3(x) + (-0.4632814174552221) * l5(x) + (-0.5506359890315792) * l7(x) + (0.6845250414579537) * l9(x) )
#define eigen_laplace_p10_6(x)   ((0.0398632115160597) * l2(x) + (0.3505633322435865) * l4(x) + (0.6448015672966024) * l6(x) + (-0.6440743612981517) * l8(x) + (0.2119326081189992) * l10(x) )
#define eigen_laplace_p10_7(x)   ((-0.2405135991654618) * l3(x) + (-0.7896100977856557) * l5(x) + (0.5482950511284399) * l7(x) + (-0.1343191683992436) * l9(x) )
#define eigen_laplace_p10_8(x)   ((0.1103045847711840) * l2(x) + (0.8971454717168891) * l4(x) + (-0.4217942554423967) * l6(x) + (0.0708050561509859) * l8(x) + (-0.0062570964181706) * l10(x) )
#define eigen_laplace_p10_9(x)   ((0.9612171512437470) * l3(x) + (-0.2745680046649570) * l5(x) + (0.0259321857669735) * l7(x) + (-0.0012331709159256) * l9(x) )
#define eigen_laplace_p10_10(x)   ((-0.9927408002342284) * l2(x) + (0.1202094754863907) * l4(x) + (-0.0039219943017710) * l6(x) + (0.0000592791732471) * l8(x) )

// derivatives

// p = 2
#define deigen_laplace_p2_2(x)   ((1.0000000000000000) * dl2(x) )


// p  3
#define deigen_laplace_p3_2(x)   ((1.0000000000000000) * dl3(x) )
#define deigen_laplace_p3_3(x)   ((1.0000000000000000) * dl2(x) )


// p = 4
#define deigen_laplace_p4_2(x)   ((0.1200767132565378) * dl2(x) + (0.9927646160764934) * dl4(x) )
#define deigen_laplace_p4_3(x)   ((-0.9999999999999999) * dl3(x) )
#define deigen_laplace_p4_4(x)   ((0.9927646160764932) * dl2(x) + (-0.1200767132565379) * dl4(x) )


// p = 5
#define deigen_laplace_p5_2(x)   ((-0.2721806201771906) * dl3(x) + (-0.9622461795195447) * dl5(x) )
#define deigen_laplace_p5_3(x)   ((-0.1200767132565378) * dl2(x) + (-0.9927646160764934) * dl4(x) )
#define deigen_laplace_p5_4(x)   ((-0.9622461795195447) * dl3(x) + (0.2721806201771906) * dl5(x) )
#define deigen_laplace_p5_5(x)   ((-0.9927646160764934) * dl2(x) + (0.1200767132565378) * dl4(x) )


// p = 6
#define deigen_laplace_p6_2(x)   ((0.0460292214008892) * dl2(x) + (0.4098443624284846) * dl4(x) + (0.9109933640608021) * dl6(x) )
#define deigen_laplace_p6_3(x)   ((0.2721806201771907) * dl3(x) + (0.9622461795195447) * dl5(x) )
#define deigen_laplace_p6_4(x)   ((0.1111170386828487) * dl2(x) + (0.9041998053176488) * dl4(x) + (-0.4124023712078773) * dl6(x) )
#define deigen_laplace_p6_5(x)   ((0.9622461795195448) * dl3(x) + (-0.2721806201771906) * dl5(x) )
#define deigen_laplace_p6_6(x)   ((0.9927408093211348) * dl2(x) + (-0.1202094449247417) * dl4(x) + (0.0039210788443064) * dl6(x) )


// p = 7
#define deigen_laplace_p7_2(x)   ((0.1244273041355310) * dl3(x) + (0.5154918741058946) * dl5(x) + (0.8478124637656356) * dl7(x) )
#define deigen_laplace_p7_3(x)   ((-0.0460292214008892) * dl2(x) + (-0.4098443624284848) * dl4(x) + (-0.9109933640608022) * dl6(x) )
#define deigen_laplace_p7_4(x)   ((0.2461133039997421) * dl3(x) + (0.8117163656110983) * dl5(x) + (-0.5296647839845884) * dl7(x) )
#define deigen_laplace_p7_5(x)   ((0.1111170386828487) * dl2(x) + (0.9041998053176490) * dl4(x) + (-0.4124023712078774) * dl6(x) )
#define deigen_laplace_p7_6(x)   ((0.9612211439517415) * dl3(x) + (-0.2745626877962528) * dl5(x) + (0.0258697292255408) * dl7(x) )
#define deigen_laplace_p7_7(x)   ((-0.9927408093211348) * dl2(x) + (0.1202094449247416) * dl4(x) + (-0.0039210788443064) * dl6(x) )


// p = 8
#define deigen_laplace_p8_2(x)   ((0.0233791285597233) * dl2(x) + (0.2118352349060833) * dl4(x) + (0.5868165705262892) * dl6(x) + (0.7811693556174851) * dl8(x) )
#define deigen_laplace_p8_3(x)   ((0.1244273041355311) * dl3(x) + (0.5154918741058940) * dl5(x) + (0.8478124637656346) * dl7(x) )
#define deigen_laplace_p8_4(x)   ((0.0418264255244555) * dl2(x) + (0.3682779654706578) * dl4(x) + (0.6912121989219027) * dl6(x) + (-0.6203608517130881) * dl8(x) )
#define deigen_laplace_p8_5(x)   ((0.2461133039997420) * dl3(x) + (0.8117163656110977) * dl5(x) + (-0.5296647839845874) * dl7(x) )
#define deigen_laplace_p8_6(x)   ((-0.1103162273830400) * dl2(x) + (-0.8972440333626744) * dl4(x) + (0.4217305141658119) * dl6(x) + (-0.0701915094045806) * dl8(x) )
#define deigen_laplace_p8_7(x)   ((-0.9612211439517412) * dl3(x) + (0.2745626877962528) * dl5(x) + (-0.0258697292255406) * dl7(x) )
#define deigen_laplace_p8_8(x)   ((-0.9927408002351916) * dl2(x) + (0.1202094754840794) * dl4(x) + (-0.0039219942315734) * dl6(x) + (0.0000592746276296) * dl8(x) )


// p  9
#define deigen_laplace_p9_2(x)   ((-0.0679631762828228) * dl3(x) + (-0.2941064986189726) * dl5(x) + (-0.6288879599997316) * dl7(x) + (-0.7165070187423248) * dl9(x) )
#define deigen_laplace_p9_3(x)   ((-0.0233791285597233) * dl2(x) + (-0.2118352349060834) * dl4(x) + (-0.5868165705262897) * dl6(x) + (-0.7811693556174855) * dl8(x) )
#define deigen_laplace_p9_4(x)   ((-0.1166010439099156) * dl3(x) + (-0.4632814174552214) * dl5(x) + (-0.5506359890315790) * dl7(x) + (0.6845250414579539) * dl9(x) )
#define deigen_laplace_p9_5(x)   ((0.0418264255244555) * dl2(x) + (0.3682779654706578) * dl4(x) + (0.6912121989219027) * dl6(x) + (-0.6203608517130882) * dl8(x) )
#define deigen_laplace_p9_6(x)   ((0.2405135991654621) * dl3(x) + (0.7896100977856557) * dl5(x) + (-0.5482950511284391) * dl7(x) + (0.1343191683992437) * dl9(x) )
#define deigen_laplace_p9_7(x)   ((0.1103162273830400) * dl2(x) + (0.8972440333626747) * dl4(x) + (-0.4217305141658118) * dl6(x) + (0.0701915094045807) * dl8(x) )
#define deigen_laplace_p9_8(x)   ((0.9612171512437465) * dl3(x) + (-0.2745680046649571) * dl5(x) + (0.0259321857669741) * dl7(x) + (-0.0012331709159256) * dl9(x) )
#define deigen_laplace_p9_9(x)   ((-0.9927408002351915) * dl2(x) + (0.1202094754840794) * dl4(x) + (-0.0039219942315734) * dl6(x) + (0.0000592746276296) * dl8(x) )


// p = 10
#define deigen_laplace_p10_2(x)   ((-0.0135012819178048) * dl2(x) + (-0.1230733557674822) * dl4(x) + (-0.3644642014036448) * dl6(x) + (-0.6487214117012139) * dl8(x) + (-0.6565036484150675) * dl10(x) )
#define deigen_laplace_p10_3(x)   ((-0.0679631762828231) * dl3(x) + (-0.2941064986189725) * dl5(x) + (-0.6288879599997312) * dl7(x) + (-0.7165070187423248) * dl9(x) )
#define deigen_laplace_p10_4(x)   ((-0.0229617482557899) * dl2(x) + (-0.2064896698810138) * dl4(x) + (-0.5229445560893998) * dl6(x) + (-0.3991433171533719) * dl8(x) + (0.7239118578789349) * dl10(x) )
#define deigen_laplace_p10_5(x)   ((-0.1166010439099159) * dl3(x) + (-0.4632814174552221) * dl5(x) + (-0.5506359890315792) * dl7(x) + (0.6845250414579537) * dl9(x) )
#define deigen_laplace_p10_6(x)   ((0.0398632115160597) * dl2(x) + (0.3505633322435865) * dl4(x) + (0.6448015672966024) * dl6(x) + (-0.6440743612981517) * dl8(x) + (0.2119326081189992) * dl10(x) )
#define deigen_laplace_p10_7(x)   ((-0.2405135991654618) * dl3(x) + (-0.7896100977856557) * dl5(x) + (0.5482950511284399) * dl7(x) + (-0.1343191683992436) * dl9(x) )
#define deigen_laplace_p10_8(x)   ((0.1103045847711840) * dl2(x) + (0.8971454717168891) * dl4(x) + (-0.4217942554423967) * dl6(x) + (0.0708050561509859) * dl8(x) + (-0.0062570964181706) * dl10(x) )
#define deigen_laplace_p10_9(x)   ((0.9612171512437470) * dl3(x) + (-0.2745680046649570) * dl5(x) + (0.0259321857669735) * dl7(x) + (-0.0012331709159256) * dl9(x) )
#define deigen_laplace_p10_10(x)   ((-0.9927408002342284) * dl2(x) + (0.1202094754863907) * dl4(x) + (-0.0039219943017710) * dl6(x) + (0.0000592791732471) * dl8(x) )



// p = 2
#define d2eigen_laplace_p2_2(x)   ((1.0000000000000000) * d2l2(x) )


// p = 3
#define d2eigen_laplace_p3_2(x)   ((1.0000000000000000) * d2l3(x) )
#define d2eigen_laplace_p3_3(x)   ((1.0000000000000000) * d2l2(x) )


// p = 4
#define d2eigen_laplace_p4_2(x)   ((0.1200767132565378) * d2l2(x) + (0.9927646160764934) * d2l4(x) )
#define d2eigen_laplace_p4_3(x)   ((-0.9999999999999999) * d2l3(x) )
#define d2eigen_laplace_p4_4(x)   ((0.9927646160764932) * d2l2(x) + (-0.1200767132565379) * d2l4(x) )


// p = 5
#define d2eigen_laplace_p5_2(x)   ((-0.2721806201771906) * d2l3(x) + (-0.9622461795195447) * d2l5(x) )
#define d2eigen_laplace_p5_3(x)   ((-0.1200767132565378) * d2l2(x) + (-0.9927646160764934) * d2l4(x) )
#define d2eigen_laplace_p5_4(x)   ((-0.9622461795195447) * d2l3(x) + (0.2721806201771906) * d2l5(x) )
#define d2eigen_laplace_p5_5(x)   ((-0.9927646160764934) * d2l2(x) + (0.1200767132565378) * d2l4(x) )


// p = 6
#define d2eigen_laplace_p6_2(x)   ((0.0460292214008892) * d2l2(x) + (0.4098443624284846) * d2l4(x) + (0.9109933640608021) * d2l6(x) )
#define d2eigen_laplace_p6_3(x)   ((0.2721806201771907) * d2l3(x) + (0.9622461795195447) * d2l5(x) )
#define d2eigen_laplace_p6_4(x)   ((0.1111170386828487) * d2l2(x) + (0.9041998053176488) * d2l4(x) + (-0.4124023712078773) * d2l6(x) )
#define d2eigen_laplace_p6_5(x)   ((0.9622461795195448) * d2l3(x) + (-0.2721806201771906) * d2l5(x) )
#define d2eigen_laplace_p6_6(x)   ((0.9927408093211348) * d2l2(x) + (-0.1202094449247417) * d2l4(x) + (0.0039210788443064) * d2l6(x) )


// p = 7
#define d2eigen_laplace_p7_2(x)   ((0.1244273041355310) * d2l3(x) + (0.5154918741058946) * d2l5(x) + (0.8478124637656356) * d2l7(x) )
#define d2eigen_laplace_p7_3(x)   ((-0.0460292214008892) * d2l2(x) + (-0.4098443624284848) * d2l4(x) + (-0.9109933640608022) * d2l6(x) )
#define d2eigen_laplace_p7_4(x)   ((0.2461133039997421) * d2l3(x) + (0.8117163656110983) * d2l5(x) + (-0.5296647839845884) * d2l7(x) )
#define d2eigen_laplace_p7_5(x)   ((0.1111170386828487) * d2l2(x) + (0.9041998053176490) * d2l4(x) + (-0.4124023712078774) * d2l6(x) )
#define d2eigen_laplace_p7_6(x)   ((0.9612211439517415) * d2l3(x) + (-0.2745626877962528) * d2l5(x) + (0.0258697292255408) * d2l7(x) )
#define d2eigen_laplace_p7_7(x)   ((-0.9927408093211348) * d2l2(x) + (0.1202094449247416) * d2l4(x) + (-0.0039210788443064) * d2l6(x) )


// p = 8
#define d2eigen_laplace_p8_2(x)   ((0.0233791285597233) * d2l2(x) + (0.2118352349060833) * d2l4(x) + (0.5868165705262892) * d2l6(x) + (0.7811693556174851) * d2l8(x) )
#define d2eigen_laplace_p8_3(x)   ((0.1244273041355311) * d2l3(x) + (0.5154918741058940) * d2l5(x) + (0.8478124637656346) * d2l7(x) )
#define d2eigen_laplace_p8_4(x)   ((0.0418264255244555) * d2l2(x) + (0.3682779654706578) * d2l4(x) + (0.6912121989219027) * d2l6(x) + (-0.6203608517130881) * d2l8(x) )
#define d2eigen_laplace_p8_5(x)   ((0.2461133039997420) * d2l3(x) + (0.8117163656110977) * d2l5(x) + (-0.5296647839845874) * d2l7(x) )
#define d2eigen_laplace_p8_6(x)   ((-0.1103162273830400) * d2l2(x) + (-0.8972440333626744) * d2l4(x) + (0.4217305141658119) * d2l6(x) + (-0.0701915094045806) * d2l8(x) )
#define d2eigen_laplace_p8_7(x)   ((-0.9612211439517412) * d2l3(x) + (0.2745626877962528) * d2l5(x) + (-0.0258697292255406) * d2l7(x) )
#define d2eigen_laplace_p8_8(x)   ((-0.9927408002351916) * d2l2(x) + (0.1202094754840794) * d2l4(x) + (-0.0039219942315734) * d2l6(x) + (0.0000592746276296) * d2l8(x) )


// p = 9
#define d2eigen_laplace_p9_2(x)   ((-0.0679631762828228) * d2l3(x) + (-0.2941064986189726) * d2l5(x) + (-0.6288879599997316) * d2l7(x) + (-0.7165070187423248) * d2l9(x) )
#define d2eigen_laplace_p9_3(x)   ((-0.0233791285597233) * d2l2(x) + (-0.2118352349060834) * d2l4(x) + (-0.5868165705262897) * d2l6(x) + (-0.7811693556174855) * d2l8(x) )
#define d2eigen_laplace_p9_4(x)   ((-0.1166010439099156) * d2l3(x) + (-0.4632814174552214) * d2l5(x) + (-0.5506359890315790) * d2l7(x) + (0.6845250414579539) * d2l9(x) )
#define d2eigen_laplace_p9_5(x)   ((0.0418264255244555) * d2l2(x) + (0.3682779654706578) * d2l4(x) + (0.6912121989219027) * d2l6(x) + (-0.6203608517130882) * d2l8(x) )
#define d2eigen_laplace_p9_6(x)   ((0.2405135991654621) * d2l3(x) + (0.7896100977856557) * d2l5(x) + (-0.5482950511284391) * d2l7(x) + (0.1343191683992437) * d2l9(x) )
#define d2eigen_laplace_p9_7(x)   ((0.1103162273830400) * d2l2(x) + (0.8972440333626747) * d2l4(x) + (-0.4217305141658118) * d2l6(x) + (0.0701915094045807) * d2l8(x) )
#define d2eigen_laplace_p9_8(x)   ((0.9612171512437465) * d2l3(x) + (-0.2745680046649571) * d2l5(x) + (0.0259321857669741) * d2l7(x) + (-0.0012331709159256) * d2l9(x) )
#define d2eigen_laplace_p9_9(x)   ((-0.9927408002351915) * d2l2(x) + (0.1202094754840794) * d2l4(x) + (-0.0039219942315734) * d2l6(x) + (0.0000592746276296) * d2l8(x) )


// p = 10
#define d2eigen_laplace_p10_2(x)   ((-0.0135012819178048) * d2l2(x) + (-0.1230733557674822) * d2l4(x) + (-0.3644642014036448) * d2l6(x) + (-0.6487214117012139) * d2l8(x) + (-0.6565036484150675) * d2l10(x) )
#define d2eigen_laplace_p10_3(x)   ((-0.0679631762828231) * d2l3(x) + (-0.2941064986189725) * d2l5(x) + (-0.6288879599997312) * d2l7(x) + (-0.7165070187423248) * d2l9(x) )
#define d2eigen_laplace_p10_4(x)   ((-0.0229617482557899) * d2l2(x) + (-0.2064896698810138) * d2l4(x) + (-0.5229445560893998) * d2l6(x) + (-0.3991433171533719) * d2l8(x) + (0.7239118578789349) * d2l10(x) )
#define d2eigen_laplace_p10_5(x)   ((-0.1166010439099159) * d2l3(x) + (-0.4632814174552221) * d2l5(x) + (-0.5506359890315792) * d2l7(x) + (0.6845250414579537) * d2l9(x) )
#define d2eigen_laplace_p10_6(x)   ((0.0398632115160597) * d2l2(x) + (0.3505633322435865) * d2l4(x) + (0.6448015672966024) * d2l6(x) + (-0.6440743612981517) * d2l8(x) + (0.2119326081189992) * d2l10(x) )
#define d2eigen_laplace_p10_7(x)   ((-0.2405135991654618) * d2l3(x) + (-0.7896100977856557) * d2l5(x) + (0.5482950511284399) * d2l7(x) + (-0.1343191683992436) * d2l9(x) )
#define d2eigen_laplace_p10_8(x)   ((0.1103045847711840) * d2l2(x) + (0.8971454717168891) * d2l4(x) + (-0.4217942554423967) * d2l6(x) + (0.0708050561509859) * d2l8(x) + (-0.0062570964181706) * d2l10(x) )
#define d2eigen_laplace_p10_9(x)   ((0.9612171512437470) * d2l3(x) + (-0.2745680046649570) * d2l5(x) + (0.0259321857669735) * d2l7(x) + (-0.0012331709159256) * d2l9(x) )
#define d2eigen_laplace_p10_10(x)   ((-0.9927408002342284) * d2l2(x) + (0.1202094754863907) * d2l4(x) + (-0.0039219943017710) * d2l6(x) + (0.0000592791732471) * d2l8(x) )



////////////////////////////////////////////////////////////////////////////////////////////////////

//// quad eigen lobatto shapeset /////////////////////////////////////////////////////////////////


static double eigen_quad_0_0(double x, double y)
{
  return  l0(x) * l0(y);
}

static double eigen_quad_0_0x(double x, double y)
{
  return  dl0(x) * l0(y);
}

static double eigen_quad_0_0y(double x, double y)
{
  return  l0(x) * dl0(y);
}

static double eigen_quad_0_1(double x, double y)
{
  return  l0(x) * l1(y);
}

static double eigen_quad_0_1x(double x, double y)
{
  return  dl0(x) * l1(y);
}

static double eigen_quad_0_1y(double x, double y)
{
  return  l0(x) * dl1(y);
}

static double eigen_quad_0_2(double x, double y)
{
  return  l0(x) * l2(y);
}

static double eigen_quad_0_2x(double x, double y)
{
  return  dl0(x) * l2(y);
}

static double eigen_quad_0_2y(double x, double y)
{
  return  l0(x) * dl2(y);
}

static double eigen_quad_0_3_0(double x, double y)
{
  return -l0(x) * l3(y);
}

static double eigen_quad_0_3_1(double x, double y)
{
  return -(-l0(x) * l3(y));
}

static double eigen_quad_0_3x_0(double x, double y)
{
  return -dl0(x) * l3(y);
}

static double eigen_quad_0_3y_0(double x, double y)
{
  return -l0(x) * dl3(y);
}

static double eigen_quad_0_3x_1(double x, double y)
{
  return -(-dl0(x) * l3(y));
}

static double eigen_quad_0_3y_1(double x, double y)
{
  return -(-l0(x) * dl3(y));
}

static double eigen_quad_0_4(double x, double y)
{
  return  l0(x) * l4(y);
}

static double eigen_quad_0_4x(double x, double y)
{
  return  dl0(x) * l4(y);
}

static double eigen_quad_0_4y(double x, double y)
{
  return  l0(x) * dl4(y);
}

static double eigen_quad_0_5_0(double x, double y)
{
  return -l0(x) * l5(y);
}

static double eigen_quad_0_5_1(double x, double y)
{
  return -(-l0(x) * l5(y));
}

static double eigen_quad_0_5x_0(double x, double y)
{
  return -dl0(x) * l5(y);
}

static double eigen_quad_0_5y_0(double x, double y)
{
  return -l0(x) * dl5(y);
}

static double eigen_quad_0_5x_1(double x, double y)
{
  return -(-dl0(x) * l5(y));
}

static double eigen_quad_0_5y_1(double x, double y)
{
  return -(-l0(x) * dl5(y));
}

static double eigen_quad_0_6(double x, double y)
{
  return  l0(x) * l6(y);
}

static double eigen_quad_0_6x(double x, double y)
{
  return  dl0(x) * l6(y);
}

static double eigen_quad_0_6y(double x, double y)
{
  return  l0(x) * dl6(y);
}

static double eigen_quad_0_7_0(double x, double y)
{
  return -l0(x) * l7(y);
}

static double eigen_quad_0_7_1(double x, double y)
{
  return -(-l0(x) * l7(y));
}

static double eigen_quad_0_7x_0(double x, double y)
{
  return -dl0(x) * l7(y);
}

static double eigen_quad_0_7y_0(double x, double y)
{
  return -l0(x) * dl7(y);
}

static double eigen_quad_0_7x_1(double x, double y)
{
  return -(-dl0(x) * l7(y));
}

static double eigen_quad_0_7y_1(double x, double y)
{
  return -(-l0(x) * dl7(y));
}

static double eigen_quad_0_8(double x, double y)
{
  return  l0(x) * l8(y);
}

static double eigen_quad_0_8x(double x, double y)
{
  return  dl0(x) * l8(y);
}

static double eigen_quad_0_8y(double x, double y)
{
  return  l0(x) * dl8(y);
}

static double eigen_quad_0_9_0(double x, double y)
{
  return -l0(x) * l9(y);
}

static double eigen_quad_0_9_1(double x, double y)
{
  return -(-l0(x) * l9(y));
}

static double eigen_quad_0_9x_0(double x, double y)
{
  return -dl0(x) * l9(y);
}

static double eigen_quad_0_9y_0(double x, double y)
{
  return -l0(x) * dl9(y);
}

static double eigen_quad_0_9x_1(double x, double y)
{
  return -(-dl0(x) * l9(y));
}

static double eigen_quad_0_9y_1(double x, double y)
{
  return -(-l0(x) * dl9(y));
}

static double eigen_quad_0_10(double x, double y)
{
  return  l0(x) * l10(y);
}

static double eigen_quad_0_10x(double x, double y)
{
  return  dl0(x) * l10(y);
}

static double eigen_quad_0_10y(double x, double y)
{
  return  l0(x) * dl10(y);
}

static double eigen_quad_1_0(double x, double y)
{
  return  l1(x) * l0(y);
}

static double eigen_quad_1_0x(double x, double y)
{
  return  dl1(x) * l0(y);
}

static double eigen_quad_1_0y(double x, double y)
{
  return  l1(x) * dl0(y);
}

static double eigen_quad_1_1(double x, double y)
{
  return  l1(x) * l1(y);
}

static double eigen_quad_1_1x(double x, double y)
{
  return  dl1(x) * l1(y);
}

static double eigen_quad_1_1y(double x, double y)
{
  return  l1(x) * dl1(y);
}

static double eigen_quad_1_2(double x, double y)
{
  return  l1(x) * l2(y);
}

static double eigen_quad_1_2x(double x, double y)
{
  return  dl1(x) * l2(y);
}

static double eigen_quad_1_2y(double x, double y)
{
  return  l1(x) * dl2(y);
}

static double eigen_quad_1_3_0(double x, double y)
{
  return  l1(x) * l3(y);
}

static double eigen_quad_1_3_1(double x, double y)
{
  return -( l1(x) * l3(y));
}

static double eigen_quad_1_3x_0(double x, double y)
{
  return  dl1(x) * l3(y);
}

static double eigen_quad_1_3y_0(double x, double y)
{
  return  l1(x) * dl3(y);
}

static double eigen_quad_1_3x_1(double x, double y)
{
  return -( dl1(x) * l3(y));
}

static double eigen_quad_1_3y_1(double x, double y)
{
  return -( l1(x) * dl3(y));
}

static double eigen_quad_1_4(double x, double y)
{
  return  l1(x) * l4(y);
}

static double eigen_quad_1_4x(double x, double y)
{
  return  dl1(x) * l4(y);
}

static double eigen_quad_1_4y(double x, double y)
{
  return  l1(x) * dl4(y);
}

static double eigen_quad_1_5_0(double x, double y)
{
  return  l1(x) * l5(y);
}

static double eigen_quad_1_5_1(double x, double y)
{
  return -( l1(x) * l5(y));
}

static double eigen_quad_1_5x_0(double x, double y)
{
  return  dl1(x) * l5(y);
}

static double eigen_quad_1_5y_0(double x, double y)
{
  return  l1(x) * dl5(y);
}

static double eigen_quad_1_5x_1(double x, double y)
{
  return -( dl1(x) * l5(y));
}

static double eigen_quad_1_5y_1(double x, double y)
{
  return -( l1(x) * dl5(y));
}

static double eigen_quad_1_6(double x, double y)
{
  return  l1(x) * l6(y);
}

static double eigen_quad_1_6x(double x, double y)
{
  return  dl1(x) * l6(y);
}

static double eigen_quad_1_6y(double x, double y)
{
  return  l1(x) * dl6(y);
}

static double eigen_quad_1_7_0(double x, double y)
{
  return  l1(x) * l7(y);
}

static double eigen_quad_1_7_1(double x, double y)
{
  return -( l1(x) * l7(y));
}

static double eigen_quad_1_7x_0(double x, double y)
{
  return  dl1(x) * l7(y);
}

static double eigen_quad_1_7y_0(double x, double y)
{
  return  l1(x) * dl7(y);
}

static double eigen_quad_1_7x_1(double x, double y)
{
  return -( dl1(x) * l7(y));
}

static double eigen_quad_1_7y_1(double x, double y)
{
  return -( l1(x) * dl7(y));
}

static double eigen_quad_1_8(double x, double y)
{
  return  l1(x) * l8(y);
}

static double eigen_quad_1_8x(double x, double y)
{
  return  dl1(x) * l8(y);
}

static double eigen_quad_1_8y(double x, double y)
{
  return  l1(x) * dl8(y);
}

static double eigen_quad_1_9_0(double x, double y)
{
  return  l1(x) * l9(y);
}

static double eigen_quad_1_9_1(double x, double y)
{
  return -( l1(x) * l9(y));
}

static double eigen_quad_1_9x_0(double x, double y)
{
  return  dl1(x) * l9(y);
}

static double eigen_quad_1_9y_0(double x, double y)
{
  return  l1(x) * dl9(y);
}

static double eigen_quad_1_9x_1(double x, double y)
{
  return -( dl1(x) * l9(y));
}

static double eigen_quad_1_9y_1(double x, double y)
{
  return -( l1(x) * dl9(y));
}

static double eigen_quad_1_10(double x, double y)
{
  return  l1(x) * l10(y);
}

static double eigen_quad_1_10x(double x, double y)
{
  return  dl1(x) * l10(y);
}

static double eigen_quad_1_10y(double x, double y)
{
  return  l1(x) * dl10(y);
}

static double eigen_quad_2_0(double x, double y)
{
  return  l2(x) * l0(y);
}

static double eigen_quad_2_0x(double x, double y)
{
  return  dl2(x) * l0(y);
}

static double eigen_quad_2_0y(double x, double y)
{
  return  l2(x) * dl0(y);
}

static double eigen_quad_2_1(double x, double y)
{
  return  l2(x) * l1(y);
}

static double eigen_quad_2_1x(double x, double y)
{
  return  dl2(x) * l1(y);
}

static double eigen_quad_2_1y(double x, double y)
{
  return  l2(x) * dl1(y);
}

static double eigen_quad_3_0_0(double x, double y)
{
  return  l3(x) * l0(y);
}

static double eigen_quad_3_0_1(double x, double y)
{
  return -( l3(x) * l0(y));
}

static double eigen_quad_3_0x_0(double x, double y)
{
  return  dl3(x) * l0(y);
}

static double eigen_quad_3_0y_0(double x, double y)
{
  return  l3(x) * dl0(y);
}

static double eigen_quad_3_0x_1(double x, double y)
{
  return -( dl3(x) * l0(y));
}

static double eigen_quad_3_0y_1(double x, double y)
{
  return -( l3(x) * dl0(y));
}

static double eigen_quad_3_1_0(double x, double y)
{
  return -l3(x) * l1(y);
}

static double eigen_quad_3_1_1(double x, double y)
{
  return -(-l3(x) * l1(y));
}

static double eigen_quad_3_1x_0(double x, double y)
{
  return -dl3(x) * l1(y);
}

static double eigen_quad_3_1y_0(double x, double y)
{
  return -l3(x) * dl1(y);
}

static double eigen_quad_3_1x_1(double x, double y)
{
  return -(-dl3(x) * l1(y));
}

static double eigen_quad_3_1y_1(double x, double y)
{
  return -(-l3(x) * dl1(y));
}

static double eigen_quad_4_0(double x, double y)
{
  return  l4(x) * l0(y);
}

static double eigen_quad_4_0x(double x, double y)
{
  return  dl4(x) * l0(y);
}

static double eigen_quad_4_0y(double x, double y)
{
  return  l4(x) * dl0(y);
}

static double eigen_quad_4_1(double x, double y)
{
  return  l4(x) * l1(y);
}

static double eigen_quad_4_1x(double x, double y)
{
  return  dl4(x) * l1(y);
}

static double eigen_quad_4_1y(double x, double y)
{
  return  l4(x) * dl1(y);
}

static double eigen_quad_5_0_0(double x, double y)
{
  return  l5(x) * l0(y);
}

static double eigen_quad_5_0_1(double x, double y)
{
  return -( l5(x) * l0(y));
}

static double eigen_quad_5_0x_0(double x, double y)
{
  return  dl5(x) * l0(y);
}

static double eigen_quad_5_0y_0(double x, double y)
{
  return  l5(x) * dl0(y);
}

static double eigen_quad_5_0x_1(double x, double y)
{
  return -( dl5(x) * l0(y));
}

static double eigen_quad_5_0y_1(double x, double y)
{
  return -( l5(x) * dl0(y));
}

static double eigen_quad_5_1_0(double x, double y)
{
  return -l5(x) * l1(y);
}

static double eigen_quad_5_1_1(double x, double y)
{
  return -(-l5(x) * l1(y));
}

static double eigen_quad_5_1x_0(double x, double y)
{
  return -dl5(x) * l1(y);
}

static double eigen_quad_5_1y_0(double x, double y)
{
  return -l5(x) * dl1(y);
}

static double eigen_quad_5_1x_1(double x, double y)
{
  return -(-dl5(x) * l1(y));
}

static double eigen_quad_5_1y_1(double x, double y)
{
  return -(-l5(x) * dl1(y));
}

static double eigen_quad_6_0(double x, double y)
{
  return  l6(x) * l0(y);
}

static double eigen_quad_6_0x(double x, double y)
{
  return  dl6(x) * l0(y);
}

static double eigen_quad_6_0y(double x, double y)
{
  return  l6(x) * dl0(y);
}

static double eigen_quad_6_1(double x, double y)
{
  return  l6(x) * l1(y);
}

static double eigen_quad_6_1x(double x, double y)
{
  return  dl6(x) * l1(y);
}

static double eigen_quad_6_1y(double x, double y)
{
  return  l6(x) * dl1(y);
}

static double eigen_quad_7_0_0(double x, double y)
{
  return  l7(x) * l0(y);
}

static double eigen_quad_7_0_1(double x, double y)
{
  return -( l7(x) * l0(y));
}

static double eigen_quad_7_0x_0(double x, double y)
{
  return  dl7(x) * l0(y);
}

static double eigen_quad_7_0y_0(double x, double y)
{
  return  l7(x) * dl0(y);
}

static double eigen_quad_7_0x_1(double x, double y)
{
  return -( dl7(x) * l0(y));
}

static double eigen_quad_7_0y_1(double x, double y)
{
  return -( l7(x) * dl0(y));
}

static double eigen_quad_7_1_0(double x, double y)
{
  return -l7(x) * l1(y);
}

static double eigen_quad_7_1_1(double x, double y)
{
  return -(-l7(x) * l1(y));
}

static double eigen_quad_7_1x_0(double x, double y)
{
  return -dl7(x) * l1(y);
}

static double eigen_quad_7_1y_0(double x, double y)
{
  return -l7(x) * dl1(y);
}

static double eigen_quad_7_1x_1(double x, double y)
{
  return -(-dl7(x) * l1(y));
}

static double eigen_quad_7_1y_1(double x, double y)
{
  return -(-l7(x) * dl1(y));
}

static double eigen_quad_8_0(double x, double y)
{
  return  l8(x) * l0(y);
}

static double eigen_quad_8_0x(double x, double y)
{
  return  dl8(x) * l0(y);
}

static double eigen_quad_8_0y(double x, double y)
{
  return  l8(x) * dl0(y);
}

static double eigen_quad_8_1(double x, double y)
{
  return  l8(x) * l1(y);
}

static double eigen_quad_8_1x(double x, double y)
{
  return  dl8(x) * l1(y);
}

static double eigen_quad_8_1y(double x, double y)
{
  return  l8(x) * dl1(y);
}

static double eigen_quad_9_0_0(double x, double y)
{
  return  l9(x) * l0(y);
}

static double eigen_quad_9_0_1(double x, double y)
{
  return -( l9(x) * l0(y));
}

static double eigen_quad_9_0x_0(double x, double y)
{
  return  dl9(x) * l0(y);
}

static double eigen_quad_9_0y_0(double x, double y)
{
  return  l9(x) * dl0(y);
}

static double eigen_quad_9_0x_1(double x, double y)
{
  return -( dl9(x) * l0(y));
}

static double eigen_quad_9_0y_1(double x, double y)
{
  return -( l9(x) * dl0(y));
}

static double eigen_quad_9_1_0(double x, double y)
{
  return -l9(x) * l1(y);
}

static double eigen_quad_9_1_1(double x, double y)
{
  return -(-l9(x) * l1(y));
}

static double eigen_quad_9_1x_0(double x, double y)
{
  return -dl9(x) * l1(y);
}

static double eigen_quad_9_1y_0(double x, double y)
{
  return -l9(x) * dl1(y);
}

static double eigen_quad_9_1x_1(double x, double y)
{
  return -(-dl9(x) * l1(y));
}

static double eigen_quad_9_1y_1(double x, double y)
{
  return -(-l9(x) * dl1(y));
}

static double eigen_quad_10_0(double x, double y)
{
  return  l10(x) * l0(y);
}

static double eigen_quad_10_0x(double x, double y)
{
  return  dl10(x) * l0(y);
}

static double eigen_quad_10_0y(double x, double y)
{
  return  l10(x) * dl0(y);
}

static double eigen_quad_10_1(double x, double y)
{
  return  l10(x) * l1(y);
}

static double eigen_quad_10_1x(double x, double y)
{
  return  dl10(x) * l1(y);
}

static double eigen_quad_10_1y(double x, double y)
{
  return  l10(x) * dl1(y);
}

/* bubble */

/* order: 2, 2 */
static double eigen_quad_p22_2_2(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p22_2_2x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p22_2_2y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p2_2(y);
}

/* order: 2, 3 */
static double eigen_quad_p23_2_2(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p23_2_2x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p23_2_2y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p23_2_3(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p23_2_3x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p23_2_3y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p3_3(y);
}

/* order: 2, 4 */
static double eigen_quad_p24_2_2(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p24_2_2x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p24_2_2y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p24_2_3(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p24_2_3x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p24_2_3y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p24_2_4(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p24_2_4x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p24_2_4y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p4_4(y);
}

/* order: 2, 5 */
static double eigen_quad_p25_2_2(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p25_2_2x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p25_2_2y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p25_2_3(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p25_2_3x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p25_2_3y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p25_2_4(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p25_2_4x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p25_2_4y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p25_2_5(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p25_2_5x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p25_2_5y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p5_5(y);
}

/* order: 2, 6 */
static double eigen_quad_p26_2_2(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p26_2_2x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p26_2_2y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p26_2_3(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p26_2_3x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p26_2_3y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p26_2_4(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p26_2_4x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p26_2_4y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p26_2_5(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p26_2_5x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p26_2_5y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p26_2_6(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p26_2_6x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p26_2_6y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p6_6(y);
}

/* order: 2, 7 */
static double eigen_quad_p27_2_2(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p27_2_2x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p27_2_2y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p27_2_3(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p27_2_3x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p27_2_3y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p27_2_4(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p27_2_4x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p27_2_4y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p27_2_5(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p27_2_5x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p27_2_5y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p27_2_6(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p27_2_6x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p27_2_6y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p27_2_7(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p27_2_7x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p27_2_7y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p7_7(y);
}

/* order: 2, 8 */
static double eigen_quad_p28_2_2(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p28_2_2x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p28_2_2y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p28_2_3(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p28_2_3x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p28_2_3y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p28_2_4(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p28_2_4x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p28_2_4y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p28_2_5(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p28_2_5x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p28_2_5y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p28_2_6(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p28_2_6x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p28_2_6y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p28_2_7(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p28_2_7x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p28_2_7y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p28_2_8(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p28_2_8x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p28_2_8y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p8_8(y);
}

/* order: 2, 9 */
static double eigen_quad_p29_2_2(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p29_2_2x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p29_2_2y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p29_2_3(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p29_2_3x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p29_2_3y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p29_2_4(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p29_2_4x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p29_2_4y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p29_2_5(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p29_2_5x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p29_2_5y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p29_2_6(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p29_2_6x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p29_2_6y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p29_2_7(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p29_2_7x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p29_2_7y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p29_2_8(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p29_2_8x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p29_2_8y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p29_2_9(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p29_2_9x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p29_2_9y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p9_9(y);
}

/* order: 2, 10 */
static double eigen_quad_p210_2_2(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p210_2_2x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p210_2_2y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p210_2_3(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p210_2_3x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p210_2_3y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p210_2_4(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p210_2_4x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p210_2_4y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p210_2_5(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p210_2_5x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p210_2_5y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p210_2_6(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p210_2_6x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p210_2_6y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p210_2_7(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p210_2_7x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p210_2_7y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p210_2_8(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p210_2_8x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p210_2_8y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p210_2_9(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p210_2_9x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p210_2_9y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p210_2_10(double x, double y)
{
  return  eigen_laplace_p2_2(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p210_2_10x(double x, double y)
{
  return  deigen_laplace_p2_2(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p210_2_10y(double x, double y)
{
  return  eigen_laplace_p2_2(x) * deigen_laplace_p10_10(y);
}

/* order: 3, 2 */
static double eigen_quad_p32_2_2(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p32_2_2x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p32_2_2y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p32_3_2(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p32_3_2x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p32_3_2y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p2_2(y);
}

/* order: 3, 3 */
static double eigen_quad_p33_2_2(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p33_2_2x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p33_2_2y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p33_2_3(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p33_2_3x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p33_2_3y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p33_3_2(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p33_3_2x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p33_3_2y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p33_3_3(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p33_3_3x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p33_3_3y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p3_3(y);
}

/* order: 3, 4 */
static double eigen_quad_p34_2_2(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p34_2_2x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p34_2_2y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p34_2_3(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p34_2_3x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p34_2_3y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p34_2_4(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p34_2_4x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p34_2_4y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p34_3_2(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p34_3_2x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p34_3_2y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p34_3_3(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p34_3_3x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p34_3_3y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p34_3_4(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p34_3_4x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p34_3_4y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p4_4(y);
}

/* order: 3, 5 */
static double eigen_quad_p35_2_2(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p35_2_2x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p35_2_2y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p35_2_3(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p35_2_3x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p35_2_3y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p35_2_4(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p35_2_4x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p35_2_4y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p35_2_5(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p35_2_5x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p35_2_5y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p35_3_2(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p35_3_2x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p35_3_2y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p35_3_3(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p35_3_3x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p35_3_3y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p35_3_4(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p35_3_4x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p35_3_4y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p35_3_5(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p35_3_5x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p35_3_5y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p5_5(y);
}

/* order: 3, 6 */
static double eigen_quad_p36_2_2(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p36_2_2x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p36_2_2y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p36_2_3(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p36_2_3x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p36_2_3y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p36_2_4(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p36_2_4x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p36_2_4y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p36_2_5(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p36_2_5x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p36_2_5y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p36_2_6(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p36_2_6x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p36_2_6y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p36_3_2(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p36_3_2x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p36_3_2y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p36_3_3(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p36_3_3x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p36_3_3y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p36_3_4(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p36_3_4x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p36_3_4y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p36_3_5(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p36_3_5x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p36_3_5y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p36_3_6(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p36_3_6x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p36_3_6y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p6_6(y);
}

/* order: 3, 7 */
static double eigen_quad_p37_2_2(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p37_2_2x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p37_2_2y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p37_2_3(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p37_2_3x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p37_2_3y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p37_2_4(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p37_2_4x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p37_2_4y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p37_2_5(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p37_2_5x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p37_2_5y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p37_2_6(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p37_2_6x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p37_2_6y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p37_2_7(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p37_2_7x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p37_2_7y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p37_3_2(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p37_3_2x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p37_3_2y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p37_3_3(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p37_3_3x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p37_3_3y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p37_3_4(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p37_3_4x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p37_3_4y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p37_3_5(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p37_3_5x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p37_3_5y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p37_3_6(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p37_3_6x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p37_3_6y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p37_3_7(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p37_3_7x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p37_3_7y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p7_7(y);
}

/* order: 3, 8 */
static double eigen_quad_p38_2_2(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p38_2_2x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p38_2_2y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p38_2_3(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p38_2_3x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p38_2_3y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p38_2_4(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p38_2_4x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p38_2_4y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p38_2_5(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p38_2_5x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p38_2_5y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p38_2_6(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p38_2_6x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p38_2_6y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p38_2_7(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p38_2_7x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p38_2_7y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p38_2_8(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p38_2_8x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p38_2_8y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p38_3_2(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p38_3_2x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p38_3_2y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p38_3_3(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p38_3_3x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p38_3_3y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p38_3_4(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p38_3_4x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p38_3_4y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p38_3_5(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p38_3_5x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p38_3_5y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p38_3_6(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p38_3_6x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p38_3_6y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p38_3_7(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p38_3_7x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p38_3_7y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p38_3_8(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p38_3_8x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p38_3_8y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p8_8(y);
}

/* order: 3, 9 */
static double eigen_quad_p39_2_2(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p39_2_2x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p39_2_2y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p39_2_3(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p39_2_3x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p39_2_3y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p39_2_4(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p39_2_4x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p39_2_4y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p39_2_5(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p39_2_5x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p39_2_5y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p39_2_6(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p39_2_6x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p39_2_6y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p39_2_7(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p39_2_7x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p39_2_7y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p39_2_8(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p39_2_8x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p39_2_8y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p39_2_9(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p39_2_9x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p39_2_9y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p39_3_2(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p39_3_2x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p39_3_2y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p39_3_3(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p39_3_3x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p39_3_3y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p39_3_4(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p39_3_4x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p39_3_4y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p39_3_5(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p39_3_5x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p39_3_5y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p39_3_6(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p39_3_6x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p39_3_6y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p39_3_7(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p39_3_7x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p39_3_7y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p39_3_8(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p39_3_8x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p39_3_8y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p39_3_9(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p39_3_9x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p39_3_9y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p9_9(y);
}

/* order: 3, 10 */
static double eigen_quad_p310_2_2(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p310_2_2x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p310_2_2y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p310_2_3(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p310_2_3x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p310_2_3y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p310_2_4(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p310_2_4x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p310_2_4y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p310_2_5(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p310_2_5x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p310_2_5y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p310_2_6(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p310_2_6x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p310_2_6y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p310_2_7(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p310_2_7x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p310_2_7y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p310_2_8(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p310_2_8x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p310_2_8y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p310_2_9(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p310_2_9x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p310_2_9y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p310_2_10(double x, double y)
{
  return  eigen_laplace_p3_2(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p310_2_10x(double x, double y)
{
  return  deigen_laplace_p3_2(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p310_2_10y(double x, double y)
{
  return  eigen_laplace_p3_2(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p310_3_2(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p310_3_2x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p310_3_2y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p310_3_3(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p310_3_3x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p310_3_3y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p310_3_4(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p310_3_4x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p310_3_4y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p310_3_5(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p310_3_5x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p310_3_5y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p310_3_6(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p310_3_6x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p310_3_6y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p310_3_7(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p310_3_7x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p310_3_7y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p310_3_8(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p310_3_8x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p310_3_8y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p310_3_9(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p310_3_9x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p310_3_9y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p310_3_10(double x, double y)
{
  return  eigen_laplace_p3_3(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p310_3_10x(double x, double y)
{
  return  deigen_laplace_p3_3(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p310_3_10y(double x, double y)
{
  return  eigen_laplace_p3_3(x) * deigen_laplace_p10_10(y);
}

/* order: 4, 2 */
static double eigen_quad_p42_2_2(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p42_2_2x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p42_2_2y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p42_3_2(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p42_3_2x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p42_3_2y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p42_4_2(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p42_4_2x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p42_4_2y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p2_2(y);
}

/* order: 4, 3 */
static double eigen_quad_p43_2_2(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p43_2_2x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p43_2_2y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p43_2_3(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p43_2_3x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p43_2_3y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p43_3_2(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p43_3_2x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p43_3_2y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p43_3_3(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p43_3_3x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p43_3_3y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p43_4_2(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p43_4_2x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p43_4_2y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p43_4_3(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p43_4_3x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p43_4_3y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p3_3(y);
}

/* order: 4, 4 */
static double eigen_quad_p44_2_2(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p44_2_2x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p44_2_2y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p44_2_3(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p44_2_3x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p44_2_3y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p44_2_4(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p44_2_4x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p44_2_4y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p44_3_2(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p44_3_2x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p44_3_2y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p44_3_3(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p44_3_3x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p44_3_3y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p44_3_4(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p44_3_4x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p44_3_4y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p44_4_2(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p44_4_2x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p44_4_2y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p44_4_3(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p44_4_3x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p44_4_3y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p44_4_4(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p44_4_4x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p44_4_4y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p4_4(y);
}

/* order: 4, 5 */
static double eigen_quad_p45_2_2(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p45_2_2x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p45_2_2y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p45_2_3(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p45_2_3x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p45_2_3y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p45_2_4(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p45_2_4x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p45_2_4y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p45_2_5(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p45_2_5x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p45_2_5y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p45_3_2(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p45_3_2x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p45_3_2y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p45_3_3(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p45_3_3x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p45_3_3y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p45_3_4(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p45_3_4x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p45_3_4y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p45_3_5(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p45_3_5x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p45_3_5y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p45_4_2(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p45_4_2x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p45_4_2y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p45_4_3(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p45_4_3x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p45_4_3y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p45_4_4(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p45_4_4x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p45_4_4y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p45_4_5(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p45_4_5x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p45_4_5y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p5_5(y);
}

/* order: 4, 6 */
static double eigen_quad_p46_2_2(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p46_2_2x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p46_2_2y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p46_2_3(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p46_2_3x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p46_2_3y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p46_2_4(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p46_2_4x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p46_2_4y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p46_2_5(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p46_2_5x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p46_2_5y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p46_2_6(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p46_2_6x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p46_2_6y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p46_3_2(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p46_3_2x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p46_3_2y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p46_3_3(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p46_3_3x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p46_3_3y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p46_3_4(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p46_3_4x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p46_3_4y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p46_3_5(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p46_3_5x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p46_3_5y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p46_3_6(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p46_3_6x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p46_3_6y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p46_4_2(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p46_4_2x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p46_4_2y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p46_4_3(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p46_4_3x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p46_4_3y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p46_4_4(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p46_4_4x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p46_4_4y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p46_4_5(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p46_4_5x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p46_4_5y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p46_4_6(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p46_4_6x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p46_4_6y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p6_6(y);
}

/* order: 4, 7 */
static double eigen_quad_p47_2_2(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p47_2_2x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p47_2_2y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p47_2_3(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p47_2_3x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p47_2_3y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p47_2_4(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p47_2_4x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p47_2_4y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p47_2_5(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p47_2_5x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p47_2_5y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p47_2_6(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p47_2_6x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p47_2_6y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p47_2_7(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p47_2_7x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p47_2_7y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p47_3_2(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p47_3_2x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p47_3_2y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p47_3_3(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p47_3_3x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p47_3_3y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p47_3_4(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p47_3_4x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p47_3_4y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p47_3_5(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p47_3_5x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p47_3_5y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p47_3_6(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p47_3_6x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p47_3_6y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p47_3_7(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p47_3_7x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p47_3_7y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p47_4_2(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p47_4_2x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p47_4_2y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p47_4_3(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p47_4_3x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p47_4_3y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p47_4_4(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p47_4_4x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p47_4_4y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p47_4_5(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p47_4_5x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p47_4_5y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p47_4_6(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p47_4_6x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p47_4_6y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p47_4_7(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p47_4_7x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p47_4_7y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p7_7(y);
}

/* order: 4, 8 */
static double eigen_quad_p48_2_2(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p48_2_2x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p48_2_2y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p48_2_3(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p48_2_3x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p48_2_3y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p48_2_4(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p48_2_4x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p48_2_4y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p48_2_5(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p48_2_5x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p48_2_5y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p48_2_6(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p48_2_6x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p48_2_6y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p48_2_7(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p48_2_7x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p48_2_7y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p48_2_8(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p48_2_8x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p48_2_8y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p48_3_2(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p48_3_2x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p48_3_2y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p48_3_3(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p48_3_3x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p48_3_3y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p48_3_4(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p48_3_4x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p48_3_4y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p48_3_5(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p48_3_5x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p48_3_5y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p48_3_6(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p48_3_6x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p48_3_6y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p48_3_7(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p48_3_7x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p48_3_7y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p48_3_8(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p48_3_8x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p48_3_8y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p48_4_2(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p48_4_2x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p48_4_2y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p48_4_3(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p48_4_3x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p48_4_3y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p48_4_4(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p48_4_4x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p48_4_4y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p48_4_5(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p48_4_5x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p48_4_5y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p48_4_6(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p48_4_6x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p48_4_6y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p48_4_7(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p48_4_7x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p48_4_7y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p48_4_8(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p48_4_8x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p48_4_8y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p8_8(y);
}

/* order: 4, 9 */
static double eigen_quad_p49_2_2(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p49_2_2x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p49_2_2y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p49_2_3(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p49_2_3x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p49_2_3y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p49_2_4(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p49_2_4x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p49_2_4y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p49_2_5(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p49_2_5x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p49_2_5y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p49_2_6(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p49_2_6x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p49_2_6y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p49_2_7(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p49_2_7x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p49_2_7y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p49_2_8(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p49_2_8x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p49_2_8y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p49_2_9(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p49_2_9x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p49_2_9y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p49_3_2(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p49_3_2x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p49_3_2y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p49_3_3(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p49_3_3x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p49_3_3y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p49_3_4(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p49_3_4x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p49_3_4y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p49_3_5(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p49_3_5x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p49_3_5y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p49_3_6(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p49_3_6x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p49_3_6y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p49_3_7(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p49_3_7x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p49_3_7y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p49_3_8(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p49_3_8x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p49_3_8y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p49_3_9(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p49_3_9x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p49_3_9y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p49_4_2(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p49_4_2x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p49_4_2y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p49_4_3(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p49_4_3x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p49_4_3y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p49_4_4(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p49_4_4x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p49_4_4y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p49_4_5(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p49_4_5x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p49_4_5y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p49_4_6(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p49_4_6x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p49_4_6y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p49_4_7(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p49_4_7x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p49_4_7y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p49_4_8(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p49_4_8x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p49_4_8y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p49_4_9(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p49_4_9x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p49_4_9y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p9_9(y);
}

/* order: 4, 10 */
static double eigen_quad_p410_2_2(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p410_2_2x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p410_2_2y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p410_2_3(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p410_2_3x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p410_2_3y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p410_2_4(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p410_2_4x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p410_2_4y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p410_2_5(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p410_2_5x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p410_2_5y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p410_2_6(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p410_2_6x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p410_2_6y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p410_2_7(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p410_2_7x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p410_2_7y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p410_2_8(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p410_2_8x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p410_2_8y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p410_2_9(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p410_2_9x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p410_2_9y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p410_2_10(double x, double y)
{
  return  eigen_laplace_p4_2(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p410_2_10x(double x, double y)
{
  return  deigen_laplace_p4_2(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p410_2_10y(double x, double y)
{
  return  eigen_laplace_p4_2(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p410_3_2(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p410_3_2x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p410_3_2y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p410_3_3(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p410_3_3x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p410_3_3y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p410_3_4(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p410_3_4x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p410_3_4y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p410_3_5(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p410_3_5x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p410_3_5y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p410_3_6(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p410_3_6x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p410_3_6y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p410_3_7(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p410_3_7x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p410_3_7y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p410_3_8(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p410_3_8x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p410_3_8y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p410_3_9(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p410_3_9x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p410_3_9y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p410_3_10(double x, double y)
{
  return  eigen_laplace_p4_3(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p410_3_10x(double x, double y)
{
  return  deigen_laplace_p4_3(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p410_3_10y(double x, double y)
{
  return  eigen_laplace_p4_3(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p410_4_2(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p410_4_2x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p410_4_2y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p410_4_3(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p410_4_3x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p410_4_3y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p410_4_4(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p410_4_4x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p410_4_4y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p410_4_5(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p410_4_5x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p410_4_5y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p410_4_6(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p410_4_6x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p410_4_6y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p410_4_7(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p410_4_7x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p410_4_7y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p410_4_8(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p410_4_8x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p410_4_8y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p410_4_9(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p410_4_9x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p410_4_9y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p410_4_10(double x, double y)
{
  return  eigen_laplace_p4_4(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p410_4_10x(double x, double y)
{
  return  deigen_laplace_p4_4(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p410_4_10y(double x, double y)
{
  return  eigen_laplace_p4_4(x) * deigen_laplace_p10_10(y);
}

/* order: 5, 2 */
static double eigen_quad_p52_2_2(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p52_2_2x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p52_2_2y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p52_3_2(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p52_3_2x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p52_3_2y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p52_4_2(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p52_4_2x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p52_4_2y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p52_5_2(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p52_5_2x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p52_5_2y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p2_2(y);
}

/* order: 5, 3 */
static double eigen_quad_p53_2_2(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p53_2_2x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p53_2_2y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p53_2_3(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p53_2_3x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p53_2_3y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p53_3_2(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p53_3_2x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p53_3_2y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p53_3_3(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p53_3_3x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p53_3_3y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p53_4_2(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p53_4_2x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p53_4_2y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p53_4_3(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p53_4_3x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p53_4_3y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p53_5_2(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p53_5_2x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p53_5_2y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p53_5_3(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p53_5_3x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p53_5_3y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p3_3(y);
}

/* order: 5, 4 */
static double eigen_quad_p54_2_2(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p54_2_2x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p54_2_2y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p54_2_3(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p54_2_3x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p54_2_3y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p54_2_4(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p54_2_4x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p54_2_4y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p54_3_2(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p54_3_2x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p54_3_2y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p54_3_3(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p54_3_3x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p54_3_3y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p54_3_4(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p54_3_4x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p54_3_4y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p54_4_2(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p54_4_2x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p54_4_2y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p54_4_3(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p54_4_3x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p54_4_3y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p54_4_4(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p54_4_4x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p54_4_4y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p54_5_2(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p54_5_2x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p54_5_2y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p54_5_3(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p54_5_3x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p54_5_3y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p54_5_4(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p54_5_4x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p54_5_4y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p4_4(y);
}

/* order: 5, 5 */
static double eigen_quad_p55_2_2(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p55_2_2x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p55_2_2y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p55_2_3(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p55_2_3x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p55_2_3y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p55_2_4(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p55_2_4x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p55_2_4y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p55_2_5(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p55_2_5x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p55_2_5y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p55_3_2(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p55_3_2x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p55_3_2y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p55_3_3(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p55_3_3x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p55_3_3y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p55_3_4(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p55_3_4x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p55_3_4y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p55_3_5(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p55_3_5x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p55_3_5y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p55_4_2(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p55_4_2x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p55_4_2y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p55_4_3(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p55_4_3x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p55_4_3y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p55_4_4(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p55_4_4x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p55_4_4y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p55_4_5(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p55_4_5x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p55_4_5y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p55_5_2(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p55_5_2x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p55_5_2y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p55_5_3(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p55_5_3x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p55_5_3y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p55_5_4(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p55_5_4x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p55_5_4y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p55_5_5(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p55_5_5x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p55_5_5y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p5_5(y);
}

/* order: 5, 6 */
static double eigen_quad_p56_2_2(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p56_2_2x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p56_2_2y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p56_2_3(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p56_2_3x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p56_2_3y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p56_2_4(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p56_2_4x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p56_2_4y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p56_2_5(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p56_2_5x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p56_2_5y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p56_2_6(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p56_2_6x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p56_2_6y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p56_3_2(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p56_3_2x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p56_3_2y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p56_3_3(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p56_3_3x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p56_3_3y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p56_3_4(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p56_3_4x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p56_3_4y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p56_3_5(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p56_3_5x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p56_3_5y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p56_3_6(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p56_3_6x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p56_3_6y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p56_4_2(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p56_4_2x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p56_4_2y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p56_4_3(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p56_4_3x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p56_4_3y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p56_4_4(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p56_4_4x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p56_4_4y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p56_4_5(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p56_4_5x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p56_4_5y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p56_4_6(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p56_4_6x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p56_4_6y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p56_5_2(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p56_5_2x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p56_5_2y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p56_5_3(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p56_5_3x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p56_5_3y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p56_5_4(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p56_5_4x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p56_5_4y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p56_5_5(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p56_5_5x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p56_5_5y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p56_5_6(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p56_5_6x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p56_5_6y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p6_6(y);
}

/* order: 5, 7 */
static double eigen_quad_p57_2_2(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p57_2_2x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p57_2_2y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p57_2_3(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p57_2_3x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p57_2_3y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p57_2_4(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p57_2_4x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p57_2_4y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p57_2_5(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p57_2_5x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p57_2_5y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p57_2_6(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p57_2_6x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p57_2_6y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p57_2_7(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p57_2_7x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p57_2_7y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p57_3_2(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p57_3_2x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p57_3_2y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p57_3_3(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p57_3_3x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p57_3_3y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p57_3_4(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p57_3_4x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p57_3_4y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p57_3_5(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p57_3_5x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p57_3_5y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p57_3_6(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p57_3_6x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p57_3_6y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p57_3_7(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p57_3_7x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p57_3_7y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p57_4_2(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p57_4_2x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p57_4_2y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p57_4_3(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p57_4_3x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p57_4_3y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p57_4_4(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p57_4_4x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p57_4_4y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p57_4_5(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p57_4_5x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p57_4_5y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p57_4_6(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p57_4_6x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p57_4_6y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p57_4_7(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p57_4_7x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p57_4_7y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p57_5_2(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p57_5_2x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p57_5_2y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p57_5_3(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p57_5_3x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p57_5_3y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p57_5_4(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p57_5_4x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p57_5_4y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p57_5_5(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p57_5_5x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p57_5_5y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p57_5_6(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p57_5_6x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p57_5_6y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p57_5_7(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p57_5_7x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p57_5_7y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p7_7(y);
}

/* order: 5, 8 */
static double eigen_quad_p58_2_2(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p58_2_2x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p58_2_2y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p58_2_3(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p58_2_3x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p58_2_3y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p58_2_4(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p58_2_4x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p58_2_4y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p58_2_5(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p58_2_5x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p58_2_5y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p58_2_6(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p58_2_6x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p58_2_6y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p58_2_7(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p58_2_7x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p58_2_7y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p58_2_8(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p58_2_8x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p58_2_8y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p58_3_2(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p58_3_2x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p58_3_2y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p58_3_3(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p58_3_3x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p58_3_3y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p58_3_4(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p58_3_4x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p58_3_4y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p58_3_5(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p58_3_5x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p58_3_5y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p58_3_6(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p58_3_6x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p58_3_6y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p58_3_7(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p58_3_7x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p58_3_7y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p58_3_8(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p58_3_8x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p58_3_8y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p58_4_2(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p58_4_2x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p58_4_2y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p58_4_3(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p58_4_3x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p58_4_3y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p58_4_4(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p58_4_4x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p58_4_4y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p58_4_5(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p58_4_5x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p58_4_5y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p58_4_6(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p58_4_6x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p58_4_6y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p58_4_7(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p58_4_7x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p58_4_7y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p58_4_8(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p58_4_8x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p58_4_8y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p58_5_2(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p58_5_2x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p58_5_2y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p58_5_3(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p58_5_3x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p58_5_3y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p58_5_4(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p58_5_4x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p58_5_4y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p58_5_5(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p58_5_5x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p58_5_5y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p58_5_6(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p58_5_6x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p58_5_6y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p58_5_7(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p58_5_7x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p58_5_7y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p58_5_8(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p58_5_8x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p58_5_8y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p8_8(y);
}

/* order: 5, 9 */
static double eigen_quad_p59_2_2(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p59_2_2x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p59_2_2y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p59_2_3(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p59_2_3x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p59_2_3y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p59_2_4(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p59_2_4x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p59_2_4y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p59_2_5(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p59_2_5x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p59_2_5y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p59_2_6(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p59_2_6x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p59_2_6y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p59_2_7(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p59_2_7x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p59_2_7y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p59_2_8(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p59_2_8x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p59_2_8y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p59_2_9(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p59_2_9x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p59_2_9y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p59_3_2(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p59_3_2x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p59_3_2y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p59_3_3(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p59_3_3x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p59_3_3y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p59_3_4(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p59_3_4x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p59_3_4y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p59_3_5(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p59_3_5x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p59_3_5y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p59_3_6(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p59_3_6x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p59_3_6y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p59_3_7(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p59_3_7x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p59_3_7y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p59_3_8(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p59_3_8x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p59_3_8y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p59_3_9(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p59_3_9x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p59_3_9y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p59_4_2(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p59_4_2x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p59_4_2y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p59_4_3(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p59_4_3x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p59_4_3y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p59_4_4(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p59_4_4x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p59_4_4y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p59_4_5(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p59_4_5x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p59_4_5y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p59_4_6(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p59_4_6x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p59_4_6y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p59_4_7(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p59_4_7x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p59_4_7y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p59_4_8(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p59_4_8x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p59_4_8y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p59_4_9(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p59_4_9x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p59_4_9y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p59_5_2(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p59_5_2x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p59_5_2y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p59_5_3(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p59_5_3x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p59_5_3y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p59_5_4(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p59_5_4x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p59_5_4y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p59_5_5(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p59_5_5x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p59_5_5y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p59_5_6(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p59_5_6x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p59_5_6y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p59_5_7(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p59_5_7x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p59_5_7y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p59_5_8(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p59_5_8x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p59_5_8y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p59_5_9(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p59_5_9x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p59_5_9y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p9_9(y);
}

/* order: 5, 10 */
static double eigen_quad_p510_2_2(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p510_2_2x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p510_2_2y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p510_2_3(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p510_2_3x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p510_2_3y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p510_2_4(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p510_2_4x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p510_2_4y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p510_2_5(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p510_2_5x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p510_2_5y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p510_2_6(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p510_2_6x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p510_2_6y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p510_2_7(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p510_2_7x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p510_2_7y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p510_2_8(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p510_2_8x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p510_2_8y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p510_2_9(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p510_2_9x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p510_2_9y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p510_2_10(double x, double y)
{
  return  eigen_laplace_p5_2(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p510_2_10x(double x, double y)
{
  return  deigen_laplace_p5_2(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p510_2_10y(double x, double y)
{
  return  eigen_laplace_p5_2(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p510_3_2(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p510_3_2x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p510_3_2y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p510_3_3(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p510_3_3x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p510_3_3y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p510_3_4(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p510_3_4x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p510_3_4y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p510_3_5(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p510_3_5x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p510_3_5y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p510_3_6(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p510_3_6x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p510_3_6y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p510_3_7(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p510_3_7x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p510_3_7y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p510_3_8(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p510_3_8x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p510_3_8y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p510_3_9(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p510_3_9x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p510_3_9y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p510_3_10(double x, double y)
{
  return  eigen_laplace_p5_3(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p510_3_10x(double x, double y)
{
  return  deigen_laplace_p5_3(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p510_3_10y(double x, double y)
{
  return  eigen_laplace_p5_3(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p510_4_2(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p510_4_2x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p510_4_2y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p510_4_3(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p510_4_3x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p510_4_3y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p510_4_4(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p510_4_4x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p510_4_4y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p510_4_5(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p510_4_5x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p510_4_5y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p510_4_6(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p510_4_6x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p510_4_6y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p510_4_7(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p510_4_7x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p510_4_7y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p510_4_8(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p510_4_8x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p510_4_8y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p510_4_9(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p510_4_9x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p510_4_9y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p510_4_10(double x, double y)
{
  return  eigen_laplace_p5_4(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p510_4_10x(double x, double y)
{
  return  deigen_laplace_p5_4(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p510_4_10y(double x, double y)
{
  return  eigen_laplace_p5_4(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p510_5_2(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p510_5_2x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p510_5_2y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p510_5_3(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p510_5_3x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p510_5_3y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p510_5_4(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p510_5_4x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p510_5_4y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p510_5_5(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p510_5_5x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p510_5_5y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p510_5_6(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p510_5_6x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p510_5_6y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p510_5_7(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p510_5_7x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p510_5_7y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p510_5_8(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p510_5_8x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p510_5_8y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p510_5_9(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p510_5_9x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p510_5_9y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p510_5_10(double x, double y)
{
  return  eigen_laplace_p5_5(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p510_5_10x(double x, double y)
{
  return  deigen_laplace_p5_5(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p510_5_10y(double x, double y)
{
  return  eigen_laplace_p5_5(x) * deigen_laplace_p10_10(y);
}

/* order: 6, 2 */
static double eigen_quad_p62_2_2(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p62_2_2x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p62_2_2y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p62_3_2(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p62_3_2x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p62_3_2y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p62_4_2(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p62_4_2x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p62_4_2y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p62_5_2(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p62_5_2x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p62_5_2y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p62_6_2(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p62_6_2x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p62_6_2y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p2_2(y);
}

/* order: 6, 3 */
static double eigen_quad_p63_2_2(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p63_2_2x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p63_2_2y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p63_2_3(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p63_2_3x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p63_2_3y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p63_3_2(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p63_3_2x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p63_3_2y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p63_3_3(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p63_3_3x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p63_3_3y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p63_4_2(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p63_4_2x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p63_4_2y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p63_4_3(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p63_4_3x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p63_4_3y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p63_5_2(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p63_5_2x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p63_5_2y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p63_5_3(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p63_5_3x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p63_5_3y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p63_6_2(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p63_6_2x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p63_6_2y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p63_6_3(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p63_6_3x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p63_6_3y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p3_3(y);
}

/* order: 6, 4 */
static double eigen_quad_p64_2_2(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p64_2_2x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p64_2_2y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p64_2_3(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p64_2_3x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p64_2_3y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p64_2_4(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p64_2_4x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p64_2_4y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p64_3_2(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p64_3_2x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p64_3_2y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p64_3_3(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p64_3_3x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p64_3_3y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p64_3_4(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p64_3_4x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p64_3_4y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p64_4_2(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p64_4_2x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p64_4_2y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p64_4_3(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p64_4_3x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p64_4_3y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p64_4_4(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p64_4_4x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p64_4_4y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p64_5_2(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p64_5_2x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p64_5_2y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p64_5_3(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p64_5_3x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p64_5_3y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p64_5_4(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p64_5_4x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p64_5_4y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p64_6_2(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p64_6_2x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p64_6_2y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p64_6_3(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p64_6_3x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p64_6_3y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p64_6_4(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p64_6_4x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p64_6_4y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p4_4(y);
}

/* order: 6, 5 */
static double eigen_quad_p65_2_2(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p65_2_2x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p65_2_2y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p65_2_3(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p65_2_3x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p65_2_3y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p65_2_4(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p65_2_4x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p65_2_4y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p65_2_5(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p65_2_5x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p65_2_5y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p65_3_2(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p65_3_2x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p65_3_2y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p65_3_3(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p65_3_3x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p65_3_3y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p65_3_4(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p65_3_4x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p65_3_4y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p65_3_5(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p65_3_5x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p65_3_5y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p65_4_2(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p65_4_2x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p65_4_2y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p65_4_3(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p65_4_3x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p65_4_3y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p65_4_4(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p65_4_4x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p65_4_4y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p65_4_5(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p65_4_5x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p65_4_5y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p65_5_2(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p65_5_2x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p65_5_2y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p65_5_3(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p65_5_3x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p65_5_3y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p65_5_4(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p65_5_4x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p65_5_4y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p65_5_5(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p65_5_5x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p65_5_5y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p65_6_2(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p65_6_2x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p65_6_2y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p65_6_3(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p65_6_3x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p65_6_3y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p65_6_4(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p65_6_4x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p65_6_4y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p65_6_5(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p65_6_5x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p65_6_5y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p5_5(y);
}

/* order: 6, 6 */
static double eigen_quad_p66_2_2(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p66_2_2x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p66_2_2y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p66_2_3(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p66_2_3x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p66_2_3y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p66_2_4(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p66_2_4x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p66_2_4y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p66_2_5(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p66_2_5x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p66_2_5y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p66_2_6(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p66_2_6x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p66_2_6y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p66_3_2(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p66_3_2x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p66_3_2y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p66_3_3(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p66_3_3x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p66_3_3y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p66_3_4(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p66_3_4x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p66_3_4y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p66_3_5(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p66_3_5x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p66_3_5y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p66_3_6(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p66_3_6x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p66_3_6y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p66_4_2(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p66_4_2x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p66_4_2y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p66_4_3(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p66_4_3x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p66_4_3y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p66_4_4(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p66_4_4x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p66_4_4y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p66_4_5(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p66_4_5x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p66_4_5y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p66_4_6(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p66_4_6x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p66_4_6y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p66_5_2(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p66_5_2x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p66_5_2y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p66_5_3(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p66_5_3x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p66_5_3y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p66_5_4(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p66_5_4x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p66_5_4y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p66_5_5(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p66_5_5x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p66_5_5y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p66_5_6(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p66_5_6x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p66_5_6y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p66_6_2(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p66_6_2x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p66_6_2y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p66_6_3(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p66_6_3x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p66_6_3y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p66_6_4(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p66_6_4x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p66_6_4y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p66_6_5(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p66_6_5x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p66_6_5y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p66_6_6(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p66_6_6x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p66_6_6y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p6_6(y);
}

/* order: 6, 7 */
static double eigen_quad_p67_2_2(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p67_2_2x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p67_2_2y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p67_2_3(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p67_2_3x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p67_2_3y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p67_2_4(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p67_2_4x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p67_2_4y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p67_2_5(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p67_2_5x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p67_2_5y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p67_2_6(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p67_2_6x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p67_2_6y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p67_2_7(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p67_2_7x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p67_2_7y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p67_3_2(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p67_3_2x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p67_3_2y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p67_3_3(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p67_3_3x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p67_3_3y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p67_3_4(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p67_3_4x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p67_3_4y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p67_3_5(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p67_3_5x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p67_3_5y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p67_3_6(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p67_3_6x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p67_3_6y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p67_3_7(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p67_3_7x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p67_3_7y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p67_4_2(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p67_4_2x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p67_4_2y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p67_4_3(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p67_4_3x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p67_4_3y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p67_4_4(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p67_4_4x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p67_4_4y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p67_4_5(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p67_4_5x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p67_4_5y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p67_4_6(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p67_4_6x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p67_4_6y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p67_4_7(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p67_4_7x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p67_4_7y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p67_5_2(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p67_5_2x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p67_5_2y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p67_5_3(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p67_5_3x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p67_5_3y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p67_5_4(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p67_5_4x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p67_5_4y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p67_5_5(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p67_5_5x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p67_5_5y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p67_5_6(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p67_5_6x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p67_5_6y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p67_5_7(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p67_5_7x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p67_5_7y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p67_6_2(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p67_6_2x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p67_6_2y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p67_6_3(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p67_6_3x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p67_6_3y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p67_6_4(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p67_6_4x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p67_6_4y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p67_6_5(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p67_6_5x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p67_6_5y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p67_6_6(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p67_6_6x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p67_6_6y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p67_6_7(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p67_6_7x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p67_6_7y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p7_7(y);
}

/* order: 6, 8 */
static double eigen_quad_p68_2_2(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p68_2_2x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p68_2_2y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p68_2_3(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p68_2_3x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p68_2_3y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p68_2_4(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p68_2_4x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p68_2_4y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p68_2_5(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p68_2_5x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p68_2_5y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p68_2_6(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p68_2_6x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p68_2_6y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p68_2_7(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p68_2_7x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p68_2_7y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p68_2_8(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p68_2_8x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p68_2_8y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p68_3_2(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p68_3_2x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p68_3_2y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p68_3_3(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p68_3_3x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p68_3_3y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p68_3_4(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p68_3_4x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p68_3_4y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p68_3_5(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p68_3_5x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p68_3_5y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p68_3_6(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p68_3_6x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p68_3_6y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p68_3_7(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p68_3_7x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p68_3_7y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p68_3_8(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p68_3_8x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p68_3_8y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p68_4_2(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p68_4_2x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p68_4_2y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p68_4_3(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p68_4_3x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p68_4_3y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p68_4_4(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p68_4_4x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p68_4_4y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p68_4_5(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p68_4_5x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p68_4_5y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p68_4_6(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p68_4_6x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p68_4_6y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p68_4_7(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p68_4_7x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p68_4_7y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p68_4_8(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p68_4_8x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p68_4_8y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p68_5_2(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p68_5_2x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p68_5_2y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p68_5_3(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p68_5_3x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p68_5_3y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p68_5_4(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p68_5_4x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p68_5_4y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p68_5_5(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p68_5_5x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p68_5_5y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p68_5_6(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p68_5_6x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p68_5_6y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p68_5_7(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p68_5_7x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p68_5_7y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p68_5_8(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p68_5_8x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p68_5_8y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p68_6_2(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p68_6_2x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p68_6_2y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p68_6_3(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p68_6_3x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p68_6_3y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p68_6_4(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p68_6_4x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p68_6_4y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p68_6_5(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p68_6_5x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p68_6_5y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p68_6_6(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p68_6_6x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p68_6_6y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p68_6_7(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p68_6_7x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p68_6_7y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p68_6_8(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p68_6_8x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p68_6_8y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p8_8(y);
}

/* order: 6, 9 */
static double eigen_quad_p69_2_2(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p69_2_2x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p69_2_2y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p69_2_3(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p69_2_3x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p69_2_3y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p69_2_4(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p69_2_4x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p69_2_4y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p69_2_5(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p69_2_5x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p69_2_5y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p69_2_6(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p69_2_6x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p69_2_6y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p69_2_7(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p69_2_7x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p69_2_7y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p69_2_8(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p69_2_8x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p69_2_8y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p69_2_9(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p69_2_9x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p69_2_9y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p69_3_2(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p69_3_2x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p69_3_2y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p69_3_3(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p69_3_3x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p69_3_3y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p69_3_4(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p69_3_4x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p69_3_4y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p69_3_5(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p69_3_5x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p69_3_5y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p69_3_6(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p69_3_6x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p69_3_6y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p69_3_7(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p69_3_7x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p69_3_7y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p69_3_8(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p69_3_8x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p69_3_8y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p69_3_9(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p69_3_9x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p69_3_9y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p69_4_2(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p69_4_2x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p69_4_2y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p69_4_3(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p69_4_3x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p69_4_3y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p69_4_4(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p69_4_4x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p69_4_4y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p69_4_5(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p69_4_5x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p69_4_5y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p69_4_6(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p69_4_6x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p69_4_6y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p69_4_7(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p69_4_7x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p69_4_7y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p69_4_8(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p69_4_8x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p69_4_8y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p69_4_9(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p69_4_9x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p69_4_9y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p69_5_2(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p69_5_2x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p69_5_2y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p69_5_3(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p69_5_3x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p69_5_3y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p69_5_4(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p69_5_4x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p69_5_4y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p69_5_5(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p69_5_5x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p69_5_5y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p69_5_6(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p69_5_6x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p69_5_6y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p69_5_7(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p69_5_7x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p69_5_7y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p69_5_8(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p69_5_8x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p69_5_8y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p69_5_9(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p69_5_9x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p69_5_9y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p69_6_2(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p69_6_2x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p69_6_2y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p69_6_3(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p69_6_3x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p69_6_3y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p69_6_4(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p69_6_4x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p69_6_4y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p69_6_5(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p69_6_5x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p69_6_5y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p69_6_6(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p69_6_6x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p69_6_6y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p69_6_7(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p69_6_7x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p69_6_7y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p69_6_8(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p69_6_8x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p69_6_8y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p69_6_9(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p69_6_9x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p69_6_9y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p9_9(y);
}

/* order: 6, 10 */
static double eigen_quad_p610_2_2(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p610_2_2x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p610_2_2y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p610_2_3(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p610_2_3x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p610_2_3y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p610_2_4(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p610_2_4x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p610_2_4y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p610_2_5(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p610_2_5x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p610_2_5y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p610_2_6(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p610_2_6x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p610_2_6y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p610_2_7(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p610_2_7x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p610_2_7y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p610_2_8(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p610_2_8x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p610_2_8y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p610_2_9(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p610_2_9x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p610_2_9y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p610_2_10(double x, double y)
{
  return  eigen_laplace_p6_2(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p610_2_10x(double x, double y)
{
  return  deigen_laplace_p6_2(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p610_2_10y(double x, double y)
{
  return  eigen_laplace_p6_2(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p610_3_2(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p610_3_2x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p610_3_2y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p610_3_3(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p610_3_3x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p610_3_3y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p610_3_4(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p610_3_4x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p610_3_4y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p610_3_5(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p610_3_5x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p610_3_5y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p610_3_6(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p610_3_6x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p610_3_6y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p610_3_7(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p610_3_7x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p610_3_7y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p610_3_8(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p610_3_8x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p610_3_8y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p610_3_9(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p610_3_9x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p610_3_9y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p610_3_10(double x, double y)
{
  return  eigen_laplace_p6_3(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p610_3_10x(double x, double y)
{
  return  deigen_laplace_p6_3(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p610_3_10y(double x, double y)
{
  return  eigen_laplace_p6_3(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p610_4_2(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p610_4_2x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p610_4_2y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p610_4_3(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p610_4_3x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p610_4_3y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p610_4_4(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p610_4_4x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p610_4_4y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p610_4_5(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p610_4_5x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p610_4_5y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p610_4_6(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p610_4_6x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p610_4_6y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p610_4_7(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p610_4_7x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p610_4_7y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p610_4_8(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p610_4_8x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p610_4_8y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p610_4_9(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p610_4_9x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p610_4_9y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p610_4_10(double x, double y)
{
  return  eigen_laplace_p6_4(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p610_4_10x(double x, double y)
{
  return  deigen_laplace_p6_4(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p610_4_10y(double x, double y)
{
  return  eigen_laplace_p6_4(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p610_5_2(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p610_5_2x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p610_5_2y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p610_5_3(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p610_5_3x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p610_5_3y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p610_5_4(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p610_5_4x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p610_5_4y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p610_5_5(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p610_5_5x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p610_5_5y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p610_5_6(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p610_5_6x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p610_5_6y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p610_5_7(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p610_5_7x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p610_5_7y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p610_5_8(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p610_5_8x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p610_5_8y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p610_5_9(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p610_5_9x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p610_5_9y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p610_5_10(double x, double y)
{
  return  eigen_laplace_p6_5(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p610_5_10x(double x, double y)
{
  return  deigen_laplace_p6_5(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p610_5_10y(double x, double y)
{
  return  eigen_laplace_p6_5(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p610_6_2(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p610_6_2x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p610_6_2y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p610_6_3(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p610_6_3x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p610_6_3y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p610_6_4(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p610_6_4x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p610_6_4y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p610_6_5(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p610_6_5x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p610_6_5y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p610_6_6(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p610_6_6x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p610_6_6y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p610_6_7(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p610_6_7x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p610_6_7y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p610_6_8(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p610_6_8x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p610_6_8y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p610_6_9(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p610_6_9x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p610_6_9y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p610_6_10(double x, double y)
{
  return  eigen_laplace_p6_6(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p610_6_10x(double x, double y)
{
  return  deigen_laplace_p6_6(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p610_6_10y(double x, double y)
{
  return  eigen_laplace_p6_6(x) * deigen_laplace_p10_10(y);
}

/* order: 7, 2 */
static double eigen_quad_p72_2_2(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p72_2_2x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p72_2_2y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p72_3_2(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p72_3_2x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p72_3_2y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p72_4_2(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p72_4_2x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p72_4_2y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p72_5_2(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p72_5_2x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p72_5_2y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p72_6_2(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p72_6_2x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p72_6_2y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p72_7_2(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p72_7_2x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p72_7_2y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p2_2(y);
}

/* order: 7, 3 */
static double eigen_quad_p73_2_2(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p73_2_2x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p73_2_2y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p73_2_3(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p73_2_3x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p73_2_3y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p73_3_2(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p73_3_2x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p73_3_2y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p73_3_3(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p73_3_3x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p73_3_3y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p73_4_2(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p73_4_2x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p73_4_2y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p73_4_3(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p73_4_3x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p73_4_3y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p73_5_2(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p73_5_2x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p73_5_2y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p73_5_3(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p73_5_3x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p73_5_3y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p73_6_2(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p73_6_2x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p73_6_2y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p73_6_3(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p73_6_3x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p73_6_3y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p73_7_2(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p73_7_2x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p73_7_2y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p73_7_3(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p73_7_3x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p73_7_3y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p3_3(y);
}

/* order: 7, 4 */
static double eigen_quad_p74_2_2(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p74_2_2x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p74_2_2y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p74_2_3(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p74_2_3x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p74_2_3y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p74_2_4(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p74_2_4x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p74_2_4y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p74_3_2(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p74_3_2x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p74_3_2y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p74_3_3(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p74_3_3x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p74_3_3y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p74_3_4(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p74_3_4x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p74_3_4y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p74_4_2(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p74_4_2x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p74_4_2y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p74_4_3(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p74_4_3x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p74_4_3y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p74_4_4(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p74_4_4x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p74_4_4y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p74_5_2(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p74_5_2x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p74_5_2y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p74_5_3(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p74_5_3x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p74_5_3y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p74_5_4(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p74_5_4x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p74_5_4y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p74_6_2(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p74_6_2x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p74_6_2y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p74_6_3(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p74_6_3x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p74_6_3y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p74_6_4(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p74_6_4x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p74_6_4y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p74_7_2(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p74_7_2x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p74_7_2y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p74_7_3(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p74_7_3x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p74_7_3y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p74_7_4(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p74_7_4x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p74_7_4y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p4_4(y);
}

/* order: 7, 5 */
static double eigen_quad_p75_2_2(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p75_2_2x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p75_2_2y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p75_2_3(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p75_2_3x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p75_2_3y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p75_2_4(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p75_2_4x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p75_2_4y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p75_2_5(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p75_2_5x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p75_2_5y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p75_3_2(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p75_3_2x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p75_3_2y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p75_3_3(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p75_3_3x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p75_3_3y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p75_3_4(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p75_3_4x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p75_3_4y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p75_3_5(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p75_3_5x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p75_3_5y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p75_4_2(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p75_4_2x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p75_4_2y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p75_4_3(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p75_4_3x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p75_4_3y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p75_4_4(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p75_4_4x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p75_4_4y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p75_4_5(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p75_4_5x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p75_4_5y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p75_5_2(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p75_5_2x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p75_5_2y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p75_5_3(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p75_5_3x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p75_5_3y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p75_5_4(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p75_5_4x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p75_5_4y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p75_5_5(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p75_5_5x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p75_5_5y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p75_6_2(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p75_6_2x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p75_6_2y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p75_6_3(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p75_6_3x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p75_6_3y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p75_6_4(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p75_6_4x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p75_6_4y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p75_6_5(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p75_6_5x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p75_6_5y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p75_7_2(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p75_7_2x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p75_7_2y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p75_7_3(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p75_7_3x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p75_7_3y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p75_7_4(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p75_7_4x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p75_7_4y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p75_7_5(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p75_7_5x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p75_7_5y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p5_5(y);
}

/* order: 7, 6 */
static double eigen_quad_p76_2_2(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p76_2_2x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p76_2_2y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p76_2_3(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p76_2_3x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p76_2_3y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p76_2_4(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p76_2_4x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p76_2_4y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p76_2_5(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p76_2_5x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p76_2_5y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p76_2_6(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p76_2_6x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p76_2_6y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p76_3_2(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p76_3_2x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p76_3_2y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p76_3_3(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p76_3_3x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p76_3_3y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p76_3_4(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p76_3_4x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p76_3_4y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p76_3_5(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p76_3_5x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p76_3_5y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p76_3_6(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p76_3_6x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p76_3_6y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p76_4_2(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p76_4_2x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p76_4_2y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p76_4_3(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p76_4_3x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p76_4_3y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p76_4_4(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p76_4_4x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p76_4_4y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p76_4_5(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p76_4_5x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p76_4_5y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p76_4_6(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p76_4_6x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p76_4_6y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p76_5_2(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p76_5_2x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p76_5_2y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p76_5_3(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p76_5_3x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p76_5_3y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p76_5_4(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p76_5_4x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p76_5_4y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p76_5_5(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p76_5_5x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p76_5_5y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p76_5_6(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p76_5_6x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p76_5_6y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p76_6_2(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p76_6_2x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p76_6_2y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p76_6_3(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p76_6_3x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p76_6_3y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p76_6_4(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p76_6_4x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p76_6_4y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p76_6_5(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p76_6_5x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p76_6_5y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p76_6_6(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p76_6_6x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p76_6_6y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p76_7_2(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p76_7_2x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p76_7_2y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p76_7_3(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p76_7_3x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p76_7_3y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p76_7_4(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p76_7_4x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p76_7_4y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p76_7_5(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p76_7_5x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p76_7_5y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p76_7_6(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p76_7_6x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p76_7_6y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p6_6(y);
}

/* order: 7, 7 */
static double eigen_quad_p77_2_2(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p77_2_2x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p77_2_2y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p77_2_3(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p77_2_3x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p77_2_3y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p77_2_4(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p77_2_4x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p77_2_4y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p77_2_5(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p77_2_5x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p77_2_5y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p77_2_6(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p77_2_6x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p77_2_6y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p77_2_7(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p77_2_7x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p77_2_7y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p77_3_2(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p77_3_2x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p77_3_2y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p77_3_3(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p77_3_3x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p77_3_3y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p77_3_4(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p77_3_4x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p77_3_4y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p77_3_5(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p77_3_5x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p77_3_5y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p77_3_6(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p77_3_6x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p77_3_6y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p77_3_7(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p77_3_7x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p77_3_7y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p77_4_2(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p77_4_2x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p77_4_2y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p77_4_3(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p77_4_3x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p77_4_3y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p77_4_4(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p77_4_4x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p77_4_4y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p77_4_5(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p77_4_5x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p77_4_5y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p77_4_6(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p77_4_6x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p77_4_6y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p77_4_7(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p77_4_7x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p77_4_7y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p77_5_2(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p77_5_2x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p77_5_2y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p77_5_3(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p77_5_3x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p77_5_3y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p77_5_4(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p77_5_4x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p77_5_4y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p77_5_5(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p77_5_5x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p77_5_5y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p77_5_6(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p77_5_6x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p77_5_6y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p77_5_7(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p77_5_7x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p77_5_7y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p77_6_2(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p77_6_2x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p77_6_2y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p77_6_3(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p77_6_3x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p77_6_3y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p77_6_4(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p77_6_4x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p77_6_4y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p77_6_5(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p77_6_5x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p77_6_5y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p77_6_6(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p77_6_6x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p77_6_6y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p77_6_7(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p77_6_7x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p77_6_7y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p77_7_2(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p77_7_2x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p77_7_2y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p77_7_3(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p77_7_3x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p77_7_3y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p77_7_4(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p77_7_4x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p77_7_4y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p77_7_5(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p77_7_5x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p77_7_5y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p77_7_6(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p77_7_6x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p77_7_6y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p77_7_7(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p77_7_7x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p77_7_7y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p7_7(y);
}

/* order: 7, 8 */
static double eigen_quad_p78_2_2(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p78_2_2x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p78_2_2y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p78_2_3(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p78_2_3x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p78_2_3y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p78_2_4(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p78_2_4x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p78_2_4y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p78_2_5(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p78_2_5x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p78_2_5y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p78_2_6(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p78_2_6x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p78_2_6y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p78_2_7(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p78_2_7x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p78_2_7y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p78_2_8(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p78_2_8x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p78_2_8y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p78_3_2(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p78_3_2x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p78_3_2y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p78_3_3(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p78_3_3x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p78_3_3y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p78_3_4(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p78_3_4x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p78_3_4y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p78_3_5(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p78_3_5x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p78_3_5y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p78_3_6(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p78_3_6x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p78_3_6y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p78_3_7(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p78_3_7x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p78_3_7y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p78_3_8(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p78_3_8x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p78_3_8y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p78_4_2(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p78_4_2x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p78_4_2y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p78_4_3(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p78_4_3x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p78_4_3y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p78_4_4(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p78_4_4x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p78_4_4y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p78_4_5(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p78_4_5x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p78_4_5y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p78_4_6(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p78_4_6x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p78_4_6y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p78_4_7(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p78_4_7x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p78_4_7y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p78_4_8(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p78_4_8x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p78_4_8y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p78_5_2(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p78_5_2x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p78_5_2y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p78_5_3(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p78_5_3x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p78_5_3y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p78_5_4(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p78_5_4x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p78_5_4y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p78_5_5(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p78_5_5x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p78_5_5y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p78_5_6(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p78_5_6x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p78_5_6y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p78_5_7(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p78_5_7x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p78_5_7y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p78_5_8(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p78_5_8x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p78_5_8y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p78_6_2(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p78_6_2x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p78_6_2y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p78_6_3(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p78_6_3x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p78_6_3y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p78_6_4(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p78_6_4x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p78_6_4y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p78_6_5(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p78_6_5x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p78_6_5y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p78_6_6(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p78_6_6x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p78_6_6y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p78_6_7(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p78_6_7x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p78_6_7y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p78_6_8(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p78_6_8x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p78_6_8y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p78_7_2(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p78_7_2x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p78_7_2y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p78_7_3(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p78_7_3x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p78_7_3y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p78_7_4(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p78_7_4x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p78_7_4y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p78_7_5(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p78_7_5x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p78_7_5y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p78_7_6(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p78_7_6x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p78_7_6y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p78_7_7(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p78_7_7x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p78_7_7y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p78_7_8(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p78_7_8x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p78_7_8y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p8_8(y);
}

/* order: 7, 9 */
static double eigen_quad_p79_2_2(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p79_2_2x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p79_2_2y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p79_2_3(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p79_2_3x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p79_2_3y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p79_2_4(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p79_2_4x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p79_2_4y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p79_2_5(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p79_2_5x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p79_2_5y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p79_2_6(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p79_2_6x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p79_2_6y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p79_2_7(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p79_2_7x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p79_2_7y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p79_2_8(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p79_2_8x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p79_2_8y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p79_2_9(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p79_2_9x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p79_2_9y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p79_3_2(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p79_3_2x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p79_3_2y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p79_3_3(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p79_3_3x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p79_3_3y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p79_3_4(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p79_3_4x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p79_3_4y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p79_3_5(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p79_3_5x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p79_3_5y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p79_3_6(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p79_3_6x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p79_3_6y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p79_3_7(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p79_3_7x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p79_3_7y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p79_3_8(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p79_3_8x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p79_3_8y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p79_3_9(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p79_3_9x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p79_3_9y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p79_4_2(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p79_4_2x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p79_4_2y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p79_4_3(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p79_4_3x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p79_4_3y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p79_4_4(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p79_4_4x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p79_4_4y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p79_4_5(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p79_4_5x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p79_4_5y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p79_4_6(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p79_4_6x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p79_4_6y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p79_4_7(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p79_4_7x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p79_4_7y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p79_4_8(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p79_4_8x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p79_4_8y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p79_4_9(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p79_4_9x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p79_4_9y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p79_5_2(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p79_5_2x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p79_5_2y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p79_5_3(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p79_5_3x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p79_5_3y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p79_5_4(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p79_5_4x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p79_5_4y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p79_5_5(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p79_5_5x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p79_5_5y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p79_5_6(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p79_5_6x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p79_5_6y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p79_5_7(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p79_5_7x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p79_5_7y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p79_5_8(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p79_5_8x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p79_5_8y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p79_5_9(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p79_5_9x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p79_5_9y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p79_6_2(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p79_6_2x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p79_6_2y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p79_6_3(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p79_6_3x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p79_6_3y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p79_6_4(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p79_6_4x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p79_6_4y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p79_6_5(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p79_6_5x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p79_6_5y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p79_6_6(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p79_6_6x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p79_6_6y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p79_6_7(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p79_6_7x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p79_6_7y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p79_6_8(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p79_6_8x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p79_6_8y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p79_6_9(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p79_6_9x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p79_6_9y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p79_7_2(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p79_7_2x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p79_7_2y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p79_7_3(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p79_7_3x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p79_7_3y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p79_7_4(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p79_7_4x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p79_7_4y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p79_7_5(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p79_7_5x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p79_7_5y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p79_7_6(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p79_7_6x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p79_7_6y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p79_7_7(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p79_7_7x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p79_7_7y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p79_7_8(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p79_7_8x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p79_7_8y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p79_7_9(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p79_7_9x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p79_7_9y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p9_9(y);
}

/* order: 7, 10 */
static double eigen_quad_p710_2_2(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p710_2_2x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p710_2_2y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p710_2_3(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p710_2_3x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p710_2_3y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p710_2_4(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p710_2_4x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p710_2_4y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p710_2_5(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p710_2_5x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p710_2_5y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p710_2_6(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p710_2_6x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p710_2_6y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p710_2_7(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p710_2_7x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p710_2_7y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p710_2_8(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p710_2_8x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p710_2_8y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p710_2_9(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p710_2_9x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p710_2_9y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p710_2_10(double x, double y)
{
  return  eigen_laplace_p7_2(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p710_2_10x(double x, double y)
{
  return  deigen_laplace_p7_2(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p710_2_10y(double x, double y)
{
  return  eigen_laplace_p7_2(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p710_3_2(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p710_3_2x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p710_3_2y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p710_3_3(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p710_3_3x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p710_3_3y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p710_3_4(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p710_3_4x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p710_3_4y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p710_3_5(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p710_3_5x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p710_3_5y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p710_3_6(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p710_3_6x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p710_3_6y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p710_3_7(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p710_3_7x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p710_3_7y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p710_3_8(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p710_3_8x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p710_3_8y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p710_3_9(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p710_3_9x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p710_3_9y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p710_3_10(double x, double y)
{
  return  eigen_laplace_p7_3(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p710_3_10x(double x, double y)
{
  return  deigen_laplace_p7_3(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p710_3_10y(double x, double y)
{
  return  eigen_laplace_p7_3(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p710_4_2(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p710_4_2x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p710_4_2y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p710_4_3(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p710_4_3x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p710_4_3y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p710_4_4(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p710_4_4x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p710_4_4y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p710_4_5(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p710_4_5x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p710_4_5y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p710_4_6(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p710_4_6x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p710_4_6y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p710_4_7(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p710_4_7x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p710_4_7y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p710_4_8(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p710_4_8x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p710_4_8y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p710_4_9(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p710_4_9x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p710_4_9y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p710_4_10(double x, double y)
{
  return  eigen_laplace_p7_4(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p710_4_10x(double x, double y)
{
  return  deigen_laplace_p7_4(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p710_4_10y(double x, double y)
{
  return  eigen_laplace_p7_4(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p710_5_2(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p710_5_2x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p710_5_2y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p710_5_3(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p710_5_3x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p710_5_3y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p710_5_4(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p710_5_4x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p710_5_4y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p710_5_5(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p710_5_5x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p710_5_5y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p710_5_6(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p710_5_6x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p710_5_6y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p710_5_7(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p710_5_7x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p710_5_7y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p710_5_8(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p710_5_8x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p710_5_8y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p710_5_9(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p710_5_9x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p710_5_9y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p710_5_10(double x, double y)
{
  return  eigen_laplace_p7_5(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p710_5_10x(double x, double y)
{
  return  deigen_laplace_p7_5(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p710_5_10y(double x, double y)
{
  return  eigen_laplace_p7_5(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p710_6_2(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p710_6_2x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p710_6_2y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p710_6_3(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p710_6_3x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p710_6_3y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p710_6_4(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p710_6_4x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p710_6_4y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p710_6_5(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p710_6_5x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p710_6_5y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p710_6_6(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p710_6_6x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p710_6_6y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p710_6_7(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p710_6_7x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p710_6_7y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p710_6_8(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p710_6_8x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p710_6_8y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p710_6_9(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p710_6_9x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p710_6_9y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p710_6_10(double x, double y)
{
  return  eigen_laplace_p7_6(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p710_6_10x(double x, double y)
{
  return  deigen_laplace_p7_6(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p710_6_10y(double x, double y)
{
  return  eigen_laplace_p7_6(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p710_7_2(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p710_7_2x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p710_7_2y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p710_7_3(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p710_7_3x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p710_7_3y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p710_7_4(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p710_7_4x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p710_7_4y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p710_7_5(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p710_7_5x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p710_7_5y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p710_7_6(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p710_7_6x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p710_7_6y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p710_7_7(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p710_7_7x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p710_7_7y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p710_7_8(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p710_7_8x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p710_7_8y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p710_7_9(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p710_7_9x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p710_7_9y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p710_7_10(double x, double y)
{
  return  eigen_laplace_p7_7(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p710_7_10x(double x, double y)
{
  return  deigen_laplace_p7_7(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p710_7_10y(double x, double y)
{
  return  eigen_laplace_p7_7(x) * deigen_laplace_p10_10(y);
}

/* order: 8, 2 */
static double eigen_quad_p82_2_2(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p82_2_2x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p82_2_2y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p82_3_2(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p82_3_2x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p82_3_2y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p82_4_2(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p82_4_2x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p82_4_2y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p82_5_2(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p82_5_2x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p82_5_2y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p82_6_2(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p82_6_2x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p82_6_2y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p82_7_2(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p82_7_2x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p82_7_2y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p82_8_2(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p82_8_2x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p82_8_2y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p2_2(y);
}

/* order: 8, 3 */
static double eigen_quad_p83_2_2(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p83_2_2x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p83_2_2y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p83_2_3(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p83_2_3x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p83_2_3y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p83_3_2(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p83_3_2x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p83_3_2y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p83_3_3(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p83_3_3x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p83_3_3y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p83_4_2(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p83_4_2x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p83_4_2y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p83_4_3(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p83_4_3x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p83_4_3y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p83_5_2(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p83_5_2x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p83_5_2y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p83_5_3(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p83_5_3x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p83_5_3y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p83_6_2(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p83_6_2x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p83_6_2y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p83_6_3(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p83_6_3x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p83_6_3y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p83_7_2(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p83_7_2x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p83_7_2y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p83_7_3(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p83_7_3x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p83_7_3y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p83_8_2(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p83_8_2x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p83_8_2y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p83_8_3(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p83_8_3x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p83_8_3y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p3_3(y);
}

/* order: 8, 4 */
static double eigen_quad_p84_2_2(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p84_2_2x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p84_2_2y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p84_2_3(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p84_2_3x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p84_2_3y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p84_2_4(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p84_2_4x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p84_2_4y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p84_3_2(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p84_3_2x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p84_3_2y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p84_3_3(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p84_3_3x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p84_3_3y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p84_3_4(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p84_3_4x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p84_3_4y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p84_4_2(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p84_4_2x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p84_4_2y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p84_4_3(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p84_4_3x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p84_4_3y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p84_4_4(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p84_4_4x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p84_4_4y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p84_5_2(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p84_5_2x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p84_5_2y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p84_5_3(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p84_5_3x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p84_5_3y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p84_5_4(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p84_5_4x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p84_5_4y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p84_6_2(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p84_6_2x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p84_6_2y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p84_6_3(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p84_6_3x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p84_6_3y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p84_6_4(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p84_6_4x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p84_6_4y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p84_7_2(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p84_7_2x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p84_7_2y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p84_7_3(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p84_7_3x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p84_7_3y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p84_7_4(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p84_7_4x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p84_7_4y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p84_8_2(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p84_8_2x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p84_8_2y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p84_8_3(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p84_8_3x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p84_8_3y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p84_8_4(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p84_8_4x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p84_8_4y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p4_4(y);
}

/* order: 8, 5 */
static double eigen_quad_p85_2_2(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p85_2_2x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p85_2_2y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p85_2_3(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p85_2_3x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p85_2_3y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p85_2_4(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p85_2_4x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p85_2_4y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p85_2_5(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p85_2_5x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p85_2_5y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p85_3_2(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p85_3_2x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p85_3_2y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p85_3_3(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p85_3_3x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p85_3_3y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p85_3_4(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p85_3_4x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p85_3_4y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p85_3_5(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p85_3_5x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p85_3_5y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p85_4_2(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p85_4_2x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p85_4_2y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p85_4_3(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p85_4_3x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p85_4_3y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p85_4_4(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p85_4_4x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p85_4_4y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p85_4_5(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p85_4_5x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p85_4_5y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p85_5_2(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p85_5_2x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p85_5_2y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p85_5_3(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p85_5_3x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p85_5_3y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p85_5_4(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p85_5_4x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p85_5_4y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p85_5_5(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p85_5_5x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p85_5_5y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p85_6_2(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p85_6_2x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p85_6_2y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p85_6_3(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p85_6_3x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p85_6_3y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p85_6_4(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p85_6_4x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p85_6_4y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p85_6_5(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p85_6_5x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p85_6_5y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p85_7_2(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p85_7_2x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p85_7_2y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p85_7_3(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p85_7_3x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p85_7_3y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p85_7_4(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p85_7_4x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p85_7_4y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p85_7_5(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p85_7_5x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p85_7_5y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p85_8_2(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p85_8_2x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p85_8_2y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p85_8_3(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p85_8_3x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p85_8_3y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p85_8_4(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p85_8_4x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p85_8_4y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p85_8_5(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p85_8_5x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p85_8_5y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p5_5(y);
}

/* order: 8, 6 */
static double eigen_quad_p86_2_2(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p86_2_2x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p86_2_2y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p86_2_3(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p86_2_3x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p86_2_3y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p86_2_4(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p86_2_4x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p86_2_4y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p86_2_5(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p86_2_5x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p86_2_5y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p86_2_6(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p86_2_6x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p86_2_6y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p86_3_2(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p86_3_2x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p86_3_2y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p86_3_3(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p86_3_3x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p86_3_3y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p86_3_4(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p86_3_4x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p86_3_4y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p86_3_5(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p86_3_5x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p86_3_5y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p86_3_6(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p86_3_6x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p86_3_6y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p86_4_2(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p86_4_2x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p86_4_2y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p86_4_3(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p86_4_3x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p86_4_3y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p86_4_4(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p86_4_4x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p86_4_4y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p86_4_5(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p86_4_5x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p86_4_5y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p86_4_6(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p86_4_6x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p86_4_6y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p86_5_2(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p86_5_2x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p86_5_2y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p86_5_3(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p86_5_3x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p86_5_3y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p86_5_4(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p86_5_4x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p86_5_4y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p86_5_5(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p86_5_5x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p86_5_5y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p86_5_6(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p86_5_6x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p86_5_6y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p86_6_2(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p86_6_2x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p86_6_2y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p86_6_3(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p86_6_3x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p86_6_3y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p86_6_4(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p86_6_4x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p86_6_4y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p86_6_5(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p86_6_5x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p86_6_5y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p86_6_6(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p86_6_6x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p86_6_6y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p86_7_2(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p86_7_2x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p86_7_2y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p86_7_3(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p86_7_3x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p86_7_3y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p86_7_4(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p86_7_4x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p86_7_4y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p86_7_5(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p86_7_5x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p86_7_5y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p86_7_6(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p86_7_6x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p86_7_6y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p86_8_2(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p86_8_2x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p86_8_2y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p86_8_3(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p86_8_3x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p86_8_3y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p86_8_4(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p86_8_4x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p86_8_4y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p86_8_5(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p86_8_5x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p86_8_5y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p86_8_6(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p86_8_6x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p86_8_6y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p6_6(y);
}

/* order: 8, 7 */
static double eigen_quad_p87_2_2(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p87_2_2x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p87_2_2y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p87_2_3(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p87_2_3x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p87_2_3y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p87_2_4(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p87_2_4x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p87_2_4y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p87_2_5(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p87_2_5x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p87_2_5y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p87_2_6(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p87_2_6x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p87_2_6y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p87_2_7(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p87_2_7x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p87_2_7y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p87_3_2(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p87_3_2x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p87_3_2y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p87_3_3(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p87_3_3x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p87_3_3y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p87_3_4(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p87_3_4x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p87_3_4y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p87_3_5(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p87_3_5x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p87_3_5y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p87_3_6(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p87_3_6x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p87_3_6y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p87_3_7(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p87_3_7x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p87_3_7y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p87_4_2(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p87_4_2x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p87_4_2y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p87_4_3(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p87_4_3x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p87_4_3y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p87_4_4(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p87_4_4x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p87_4_4y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p87_4_5(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p87_4_5x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p87_4_5y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p87_4_6(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p87_4_6x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p87_4_6y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p87_4_7(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p87_4_7x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p87_4_7y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p87_5_2(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p87_5_2x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p87_5_2y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p87_5_3(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p87_5_3x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p87_5_3y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p87_5_4(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p87_5_4x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p87_5_4y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p87_5_5(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p87_5_5x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p87_5_5y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p87_5_6(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p87_5_6x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p87_5_6y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p87_5_7(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p87_5_7x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p87_5_7y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p87_6_2(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p87_6_2x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p87_6_2y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p87_6_3(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p87_6_3x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p87_6_3y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p87_6_4(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p87_6_4x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p87_6_4y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p87_6_5(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p87_6_5x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p87_6_5y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p87_6_6(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p87_6_6x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p87_6_6y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p87_6_7(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p87_6_7x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p87_6_7y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p87_7_2(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p87_7_2x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p87_7_2y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p87_7_3(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p87_7_3x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p87_7_3y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p87_7_4(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p87_7_4x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p87_7_4y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p87_7_5(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p87_7_5x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p87_7_5y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p87_7_6(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p87_7_6x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p87_7_6y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p87_7_7(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p87_7_7x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p87_7_7y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p87_8_2(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p87_8_2x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p87_8_2y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p87_8_3(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p87_8_3x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p87_8_3y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p87_8_4(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p87_8_4x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p87_8_4y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p87_8_5(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p87_8_5x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p87_8_5y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p87_8_6(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p87_8_6x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p87_8_6y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p87_8_7(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p87_8_7x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p87_8_7y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p7_7(y);
}

/* order: 8, 8 */
static double eigen_quad_p88_2_2(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p88_2_2x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p88_2_2y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p88_2_3(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p88_2_3x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p88_2_3y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p88_2_4(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p88_2_4x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p88_2_4y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p88_2_5(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p88_2_5x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p88_2_5y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p88_2_6(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p88_2_6x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p88_2_6y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p88_2_7(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p88_2_7x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p88_2_7y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p88_2_8(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p88_2_8x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p88_2_8y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p88_3_2(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p88_3_2x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p88_3_2y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p88_3_3(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p88_3_3x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p88_3_3y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p88_3_4(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p88_3_4x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p88_3_4y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p88_3_5(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p88_3_5x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p88_3_5y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p88_3_6(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p88_3_6x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p88_3_6y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p88_3_7(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p88_3_7x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p88_3_7y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p88_3_8(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p88_3_8x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p88_3_8y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p88_4_2(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p88_4_2x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p88_4_2y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p88_4_3(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p88_4_3x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p88_4_3y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p88_4_4(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p88_4_4x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p88_4_4y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p88_4_5(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p88_4_5x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p88_4_5y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p88_4_6(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p88_4_6x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p88_4_6y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p88_4_7(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p88_4_7x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p88_4_7y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p88_4_8(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p88_4_8x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p88_4_8y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p88_5_2(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p88_5_2x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p88_5_2y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p88_5_3(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p88_5_3x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p88_5_3y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p88_5_4(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p88_5_4x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p88_5_4y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p88_5_5(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p88_5_5x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p88_5_5y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p88_5_6(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p88_5_6x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p88_5_6y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p88_5_7(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p88_5_7x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p88_5_7y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p88_5_8(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p88_5_8x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p88_5_8y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p88_6_2(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p88_6_2x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p88_6_2y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p88_6_3(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p88_6_3x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p88_6_3y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p88_6_4(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p88_6_4x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p88_6_4y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p88_6_5(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p88_6_5x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p88_6_5y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p88_6_6(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p88_6_6x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p88_6_6y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p88_6_7(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p88_6_7x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p88_6_7y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p88_6_8(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p88_6_8x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p88_6_8y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p88_7_2(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p88_7_2x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p88_7_2y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p88_7_3(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p88_7_3x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p88_7_3y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p88_7_4(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p88_7_4x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p88_7_4y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p88_7_5(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p88_7_5x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p88_7_5y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p88_7_6(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p88_7_6x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p88_7_6y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p88_7_7(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p88_7_7x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p88_7_7y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p88_7_8(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p88_7_8x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p88_7_8y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p88_8_2(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p88_8_2x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p88_8_2y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p88_8_3(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p88_8_3x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p88_8_3y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p88_8_4(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p88_8_4x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p88_8_4y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p88_8_5(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p88_8_5x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p88_8_5y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p88_8_6(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p88_8_6x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p88_8_6y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p88_8_7(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p88_8_7x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p88_8_7y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p88_8_8(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p88_8_8x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p88_8_8y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p8_8(y);
}

/* order: 8, 9 */
static double eigen_quad_p89_2_2(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p89_2_2x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p89_2_2y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p89_2_3(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p89_2_3x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p89_2_3y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p89_2_4(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p89_2_4x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p89_2_4y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p89_2_5(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p89_2_5x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p89_2_5y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p89_2_6(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p89_2_6x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p89_2_6y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p89_2_7(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p89_2_7x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p89_2_7y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p89_2_8(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p89_2_8x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p89_2_8y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p89_2_9(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p89_2_9x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p89_2_9y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p89_3_2(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p89_3_2x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p89_3_2y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p89_3_3(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p89_3_3x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p89_3_3y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p89_3_4(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p89_3_4x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p89_3_4y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p89_3_5(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p89_3_5x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p89_3_5y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p89_3_6(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p89_3_6x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p89_3_6y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p89_3_7(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p89_3_7x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p89_3_7y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p89_3_8(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p89_3_8x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p89_3_8y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p89_3_9(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p89_3_9x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p89_3_9y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p89_4_2(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p89_4_2x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p89_4_2y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p89_4_3(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p89_4_3x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p89_4_3y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p89_4_4(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p89_4_4x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p89_4_4y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p89_4_5(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p89_4_5x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p89_4_5y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p89_4_6(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p89_4_6x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p89_4_6y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p89_4_7(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p89_4_7x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p89_4_7y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p89_4_8(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p89_4_8x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p89_4_8y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p89_4_9(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p89_4_9x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p89_4_9y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p89_5_2(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p89_5_2x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p89_5_2y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p89_5_3(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p89_5_3x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p89_5_3y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p89_5_4(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p89_5_4x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p89_5_4y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p89_5_5(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p89_5_5x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p89_5_5y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p89_5_6(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p89_5_6x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p89_5_6y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p89_5_7(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p89_5_7x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p89_5_7y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p89_5_8(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p89_5_8x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p89_5_8y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p89_5_9(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p89_5_9x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p89_5_9y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p89_6_2(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p89_6_2x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p89_6_2y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p89_6_3(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p89_6_3x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p89_6_3y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p89_6_4(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p89_6_4x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p89_6_4y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p89_6_5(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p89_6_5x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p89_6_5y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p89_6_6(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p89_6_6x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p89_6_6y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p89_6_7(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p89_6_7x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p89_6_7y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p89_6_8(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p89_6_8x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p89_6_8y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p89_6_9(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p89_6_9x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p89_6_9y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p89_7_2(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p89_7_2x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p89_7_2y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p89_7_3(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p89_7_3x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p89_7_3y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p89_7_4(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p89_7_4x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p89_7_4y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p89_7_5(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p89_7_5x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p89_7_5y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p89_7_6(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p89_7_6x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p89_7_6y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p89_7_7(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p89_7_7x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p89_7_7y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p89_7_8(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p89_7_8x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p89_7_8y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p89_7_9(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p89_7_9x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p89_7_9y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p89_8_2(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p89_8_2x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p89_8_2y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p89_8_3(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p89_8_3x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p89_8_3y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p89_8_4(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p89_8_4x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p89_8_4y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p89_8_5(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p89_8_5x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p89_8_5y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p89_8_6(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p89_8_6x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p89_8_6y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p89_8_7(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p89_8_7x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p89_8_7y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p89_8_8(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p89_8_8x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p89_8_8y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p89_8_9(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p89_8_9x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p89_8_9y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p9_9(y);
}

/* order: 8, 10 */
static double eigen_quad_p810_2_2(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p810_2_2x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p810_2_2y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p810_2_3(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p810_2_3x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p810_2_3y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p810_2_4(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p810_2_4x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p810_2_4y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p810_2_5(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p810_2_5x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p810_2_5y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p810_2_6(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p810_2_6x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p810_2_6y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p810_2_7(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p810_2_7x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p810_2_7y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p810_2_8(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p810_2_8x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p810_2_8y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p810_2_9(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p810_2_9x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p810_2_9y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p810_2_10(double x, double y)
{
  return  eigen_laplace_p8_2(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p810_2_10x(double x, double y)
{
  return  deigen_laplace_p8_2(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p810_2_10y(double x, double y)
{
  return  eigen_laplace_p8_2(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p810_3_2(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p810_3_2x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p810_3_2y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p810_3_3(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p810_3_3x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p810_3_3y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p810_3_4(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p810_3_4x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p810_3_4y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p810_3_5(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p810_3_5x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p810_3_5y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p810_3_6(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p810_3_6x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p810_3_6y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p810_3_7(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p810_3_7x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p810_3_7y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p810_3_8(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p810_3_8x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p810_3_8y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p810_3_9(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p810_3_9x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p810_3_9y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p810_3_10(double x, double y)
{
  return  eigen_laplace_p8_3(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p810_3_10x(double x, double y)
{
  return  deigen_laplace_p8_3(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p810_3_10y(double x, double y)
{
  return  eigen_laplace_p8_3(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p810_4_2(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p810_4_2x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p810_4_2y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p810_4_3(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p810_4_3x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p810_4_3y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p810_4_4(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p810_4_4x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p810_4_4y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p810_4_5(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p810_4_5x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p810_4_5y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p810_4_6(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p810_4_6x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p810_4_6y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p810_4_7(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p810_4_7x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p810_4_7y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p810_4_8(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p810_4_8x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p810_4_8y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p810_4_9(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p810_4_9x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p810_4_9y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p810_4_10(double x, double y)
{
  return  eigen_laplace_p8_4(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p810_4_10x(double x, double y)
{
  return  deigen_laplace_p8_4(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p810_4_10y(double x, double y)
{
  return  eigen_laplace_p8_4(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p810_5_2(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p810_5_2x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p810_5_2y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p810_5_3(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p810_5_3x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p810_5_3y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p810_5_4(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p810_5_4x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p810_5_4y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p810_5_5(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p810_5_5x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p810_5_5y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p810_5_6(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p810_5_6x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p810_5_6y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p810_5_7(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p810_5_7x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p810_5_7y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p810_5_8(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p810_5_8x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p810_5_8y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p810_5_9(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p810_5_9x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p810_5_9y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p810_5_10(double x, double y)
{
  return  eigen_laplace_p8_5(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p810_5_10x(double x, double y)
{
  return  deigen_laplace_p8_5(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p810_5_10y(double x, double y)
{
  return  eigen_laplace_p8_5(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p810_6_2(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p810_6_2x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p810_6_2y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p810_6_3(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p810_6_3x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p810_6_3y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p810_6_4(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p810_6_4x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p810_6_4y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p810_6_5(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p810_6_5x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p810_6_5y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p810_6_6(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p810_6_6x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p810_6_6y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p810_6_7(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p810_6_7x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p810_6_7y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p810_6_8(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p810_6_8x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p810_6_8y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p810_6_9(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p810_6_9x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p810_6_9y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p810_6_10(double x, double y)
{
  return  eigen_laplace_p8_6(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p810_6_10x(double x, double y)
{
  return  deigen_laplace_p8_6(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p810_6_10y(double x, double y)
{
  return  eigen_laplace_p8_6(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p810_7_2(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p810_7_2x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p810_7_2y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p810_7_3(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p810_7_3x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p810_7_3y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p810_7_4(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p810_7_4x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p810_7_4y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p810_7_5(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p810_7_5x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p810_7_5y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p810_7_6(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p810_7_6x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p810_7_6y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p810_7_7(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p810_7_7x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p810_7_7y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p810_7_8(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p810_7_8x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p810_7_8y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p810_7_9(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p810_7_9x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p810_7_9y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p810_7_10(double x, double y)
{
  return  eigen_laplace_p8_7(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p810_7_10x(double x, double y)
{
  return  deigen_laplace_p8_7(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p810_7_10y(double x, double y)
{
  return  eigen_laplace_p8_7(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p810_8_2(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p810_8_2x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p810_8_2y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p810_8_3(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p810_8_3x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p810_8_3y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p810_8_4(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p810_8_4x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p810_8_4y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p810_8_5(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p810_8_5x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p810_8_5y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p810_8_6(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p810_8_6x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p810_8_6y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p810_8_7(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p810_8_7x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p810_8_7y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p810_8_8(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p810_8_8x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p810_8_8y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p810_8_9(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p810_8_9x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p810_8_9y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p810_8_10(double x, double y)
{
  return  eigen_laplace_p8_8(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p810_8_10x(double x, double y)
{
  return  deigen_laplace_p8_8(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p810_8_10y(double x, double y)
{
  return  eigen_laplace_p8_8(x) * deigen_laplace_p10_10(y);
}

/* order: 9, 2 */
static double eigen_quad_p92_2_2(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p92_2_2x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p92_2_2y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p92_3_2(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p92_3_2x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p92_3_2y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p92_4_2(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p92_4_2x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p92_4_2y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p92_5_2(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p92_5_2x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p92_5_2y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p92_6_2(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p92_6_2x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p92_6_2y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p92_7_2(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p92_7_2x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p92_7_2y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p92_8_2(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p92_8_2x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p92_8_2y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p92_9_2(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p92_9_2x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p92_9_2y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p2_2(y);
}

/* order: 9, 3 */
static double eigen_quad_p93_2_2(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p93_2_2x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p93_2_2y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p93_2_3(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p93_2_3x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p93_2_3y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p93_3_2(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p93_3_2x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p93_3_2y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p93_3_3(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p93_3_3x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p93_3_3y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p93_4_2(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p93_4_2x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p93_4_2y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p93_4_3(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p93_4_3x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p93_4_3y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p93_5_2(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p93_5_2x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p93_5_2y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p93_5_3(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p93_5_3x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p93_5_3y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p93_6_2(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p93_6_2x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p93_6_2y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p93_6_3(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p93_6_3x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p93_6_3y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p93_7_2(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p93_7_2x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p93_7_2y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p93_7_3(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p93_7_3x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p93_7_3y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p93_8_2(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p93_8_2x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p93_8_2y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p93_8_3(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p93_8_3x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p93_8_3y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p93_9_2(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p93_9_2x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p93_9_2y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p93_9_3(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p93_9_3x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p93_9_3y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p3_3(y);
}

/* order: 9, 4 */
static double eigen_quad_p94_2_2(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p94_2_2x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p94_2_2y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p94_2_3(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p94_2_3x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p94_2_3y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p94_2_4(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p94_2_4x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p94_2_4y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p94_3_2(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p94_3_2x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p94_3_2y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p94_3_3(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p94_3_3x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p94_3_3y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p94_3_4(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p94_3_4x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p94_3_4y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p94_4_2(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p94_4_2x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p94_4_2y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p94_4_3(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p94_4_3x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p94_4_3y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p94_4_4(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p94_4_4x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p94_4_4y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p94_5_2(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p94_5_2x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p94_5_2y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p94_5_3(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p94_5_3x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p94_5_3y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p94_5_4(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p94_5_4x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p94_5_4y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p94_6_2(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p94_6_2x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p94_6_2y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p94_6_3(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p94_6_3x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p94_6_3y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p94_6_4(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p94_6_4x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p94_6_4y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p94_7_2(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p94_7_2x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p94_7_2y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p94_7_3(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p94_7_3x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p94_7_3y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p94_7_4(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p94_7_4x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p94_7_4y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p94_8_2(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p94_8_2x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p94_8_2y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p94_8_3(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p94_8_3x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p94_8_3y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p94_8_4(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p94_8_4x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p94_8_4y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p94_9_2(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p94_9_2x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p94_9_2y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p94_9_3(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p94_9_3x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p94_9_3y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p94_9_4(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p94_9_4x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p94_9_4y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p4_4(y);
}

/* order: 9, 5 */
static double eigen_quad_p95_2_2(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p95_2_2x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p95_2_2y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p95_2_3(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p95_2_3x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p95_2_3y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p95_2_4(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p95_2_4x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p95_2_4y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p95_2_5(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p95_2_5x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p95_2_5y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p95_3_2(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p95_3_2x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p95_3_2y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p95_3_3(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p95_3_3x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p95_3_3y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p95_3_4(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p95_3_4x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p95_3_4y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p95_3_5(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p95_3_5x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p95_3_5y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p95_4_2(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p95_4_2x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p95_4_2y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p95_4_3(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p95_4_3x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p95_4_3y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p95_4_4(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p95_4_4x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p95_4_4y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p95_4_5(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p95_4_5x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p95_4_5y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p95_5_2(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p95_5_2x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p95_5_2y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p95_5_3(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p95_5_3x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p95_5_3y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p95_5_4(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p95_5_4x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p95_5_4y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p95_5_5(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p95_5_5x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p95_5_5y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p95_6_2(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p95_6_2x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p95_6_2y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p95_6_3(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p95_6_3x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p95_6_3y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p95_6_4(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p95_6_4x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p95_6_4y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p95_6_5(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p95_6_5x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p95_6_5y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p95_7_2(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p95_7_2x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p95_7_2y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p95_7_3(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p95_7_3x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p95_7_3y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p95_7_4(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p95_7_4x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p95_7_4y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p95_7_5(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p95_7_5x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p95_7_5y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p95_8_2(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p95_8_2x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p95_8_2y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p95_8_3(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p95_8_3x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p95_8_3y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p95_8_4(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p95_8_4x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p95_8_4y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p95_8_5(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p95_8_5x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p95_8_5y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p95_9_2(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p95_9_2x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p95_9_2y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p95_9_3(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p95_9_3x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p95_9_3y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p95_9_4(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p95_9_4x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p95_9_4y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p95_9_5(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p95_9_5x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p95_9_5y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p5_5(y);
}

/* order: 9, 6 */
static double eigen_quad_p96_2_2(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p96_2_2x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p96_2_2y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p96_2_3(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p96_2_3x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p96_2_3y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p96_2_4(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p96_2_4x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p96_2_4y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p96_2_5(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p96_2_5x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p96_2_5y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p96_2_6(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p96_2_6x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p96_2_6y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p96_3_2(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p96_3_2x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p96_3_2y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p96_3_3(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p96_3_3x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p96_3_3y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p96_3_4(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p96_3_4x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p96_3_4y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p96_3_5(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p96_3_5x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p96_3_5y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p96_3_6(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p96_3_6x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p96_3_6y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p96_4_2(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p96_4_2x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p96_4_2y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p96_4_3(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p96_4_3x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p96_4_3y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p96_4_4(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p96_4_4x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p96_4_4y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p96_4_5(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p96_4_5x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p96_4_5y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p96_4_6(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p96_4_6x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p96_4_6y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p96_5_2(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p96_5_2x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p96_5_2y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p96_5_3(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p96_5_3x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p96_5_3y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p96_5_4(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p96_5_4x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p96_5_4y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p96_5_5(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p96_5_5x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p96_5_5y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p96_5_6(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p96_5_6x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p96_5_6y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p96_6_2(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p96_6_2x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p96_6_2y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p96_6_3(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p96_6_3x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p96_6_3y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p96_6_4(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p96_6_4x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p96_6_4y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p96_6_5(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p96_6_5x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p96_6_5y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p96_6_6(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p96_6_6x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p96_6_6y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p96_7_2(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p96_7_2x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p96_7_2y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p96_7_3(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p96_7_3x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p96_7_3y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p96_7_4(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p96_7_4x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p96_7_4y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p96_7_5(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p96_7_5x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p96_7_5y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p96_7_6(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p96_7_6x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p96_7_6y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p96_8_2(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p96_8_2x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p96_8_2y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p96_8_3(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p96_8_3x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p96_8_3y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p96_8_4(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p96_8_4x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p96_8_4y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p96_8_5(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p96_8_5x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p96_8_5y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p96_8_6(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p96_8_6x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p96_8_6y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p96_9_2(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p96_9_2x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p96_9_2y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p96_9_3(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p96_9_3x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p96_9_3y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p96_9_4(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p96_9_4x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p96_9_4y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p96_9_5(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p96_9_5x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p96_9_5y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p96_9_6(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p96_9_6x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p96_9_6y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p6_6(y);
}

/* order: 9, 7 */
static double eigen_quad_p97_2_2(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p97_2_2x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p97_2_2y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p97_2_3(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p97_2_3x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p97_2_3y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p97_2_4(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p97_2_4x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p97_2_4y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p97_2_5(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p97_2_5x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p97_2_5y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p97_2_6(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p97_2_6x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p97_2_6y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p97_2_7(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p97_2_7x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p97_2_7y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p97_3_2(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p97_3_2x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p97_3_2y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p97_3_3(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p97_3_3x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p97_3_3y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p97_3_4(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p97_3_4x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p97_3_4y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p97_3_5(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p97_3_5x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p97_3_5y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p97_3_6(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p97_3_6x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p97_3_6y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p97_3_7(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p97_3_7x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p97_3_7y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p97_4_2(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p97_4_2x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p97_4_2y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p97_4_3(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p97_4_3x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p97_4_3y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p97_4_4(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p97_4_4x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p97_4_4y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p97_4_5(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p97_4_5x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p97_4_5y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p97_4_6(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p97_4_6x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p97_4_6y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p97_4_7(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p97_4_7x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p97_4_7y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p97_5_2(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p97_5_2x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p97_5_2y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p97_5_3(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p97_5_3x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p97_5_3y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p97_5_4(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p97_5_4x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p97_5_4y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p97_5_5(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p97_5_5x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p97_5_5y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p97_5_6(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p97_5_6x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p97_5_6y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p97_5_7(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p97_5_7x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p97_5_7y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p97_6_2(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p97_6_2x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p97_6_2y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p97_6_3(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p97_6_3x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p97_6_3y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p97_6_4(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p97_6_4x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p97_6_4y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p97_6_5(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p97_6_5x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p97_6_5y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p97_6_6(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p97_6_6x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p97_6_6y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p97_6_7(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p97_6_7x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p97_6_7y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p97_7_2(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p97_7_2x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p97_7_2y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p97_7_3(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p97_7_3x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p97_7_3y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p97_7_4(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p97_7_4x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p97_7_4y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p97_7_5(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p97_7_5x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p97_7_5y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p97_7_6(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p97_7_6x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p97_7_6y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p97_7_7(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p97_7_7x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p97_7_7y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p97_8_2(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p97_8_2x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p97_8_2y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p97_8_3(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p97_8_3x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p97_8_3y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p97_8_4(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p97_8_4x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p97_8_4y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p97_8_5(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p97_8_5x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p97_8_5y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p97_8_6(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p97_8_6x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p97_8_6y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p97_8_7(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p97_8_7x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p97_8_7y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p97_9_2(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p97_9_2x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p97_9_2y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p97_9_3(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p97_9_3x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p97_9_3y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p97_9_4(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p97_9_4x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p97_9_4y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p97_9_5(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p97_9_5x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p97_9_5y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p97_9_6(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p97_9_6x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p97_9_6y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p97_9_7(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p97_9_7x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p97_9_7y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p7_7(y);
}

/* order: 9, 8 */
static double eigen_quad_p98_2_2(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p98_2_2x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p98_2_2y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p98_2_3(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p98_2_3x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p98_2_3y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p98_2_4(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p98_2_4x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p98_2_4y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p98_2_5(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p98_2_5x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p98_2_5y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p98_2_6(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p98_2_6x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p98_2_6y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p98_2_7(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p98_2_7x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p98_2_7y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p98_2_8(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p98_2_8x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p98_2_8y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p98_3_2(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p98_3_2x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p98_3_2y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p98_3_3(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p98_3_3x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p98_3_3y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p98_3_4(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p98_3_4x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p98_3_4y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p98_3_5(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p98_3_5x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p98_3_5y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p98_3_6(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p98_3_6x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p98_3_6y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p98_3_7(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p98_3_7x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p98_3_7y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p98_3_8(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p98_3_8x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p98_3_8y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p98_4_2(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p98_4_2x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p98_4_2y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p98_4_3(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p98_4_3x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p98_4_3y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p98_4_4(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p98_4_4x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p98_4_4y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p98_4_5(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p98_4_5x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p98_4_5y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p98_4_6(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p98_4_6x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p98_4_6y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p98_4_7(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p98_4_7x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p98_4_7y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p98_4_8(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p98_4_8x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p98_4_8y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p98_5_2(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p98_5_2x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p98_5_2y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p98_5_3(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p98_5_3x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p98_5_3y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p98_5_4(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p98_5_4x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p98_5_4y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p98_5_5(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p98_5_5x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p98_5_5y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p98_5_6(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p98_5_6x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p98_5_6y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p98_5_7(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p98_5_7x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p98_5_7y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p98_5_8(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p98_5_8x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p98_5_8y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p98_6_2(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p98_6_2x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p98_6_2y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p98_6_3(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p98_6_3x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p98_6_3y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p98_6_4(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p98_6_4x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p98_6_4y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p98_6_5(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p98_6_5x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p98_6_5y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p98_6_6(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p98_6_6x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p98_6_6y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p98_6_7(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p98_6_7x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p98_6_7y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p98_6_8(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p98_6_8x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p98_6_8y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p98_7_2(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p98_7_2x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p98_7_2y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p98_7_3(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p98_7_3x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p98_7_3y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p98_7_4(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p98_7_4x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p98_7_4y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p98_7_5(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p98_7_5x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p98_7_5y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p98_7_6(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p98_7_6x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p98_7_6y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p98_7_7(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p98_7_7x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p98_7_7y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p98_7_8(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p98_7_8x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p98_7_8y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p98_8_2(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p98_8_2x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p98_8_2y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p98_8_3(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p98_8_3x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p98_8_3y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p98_8_4(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p98_8_4x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p98_8_4y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p98_8_5(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p98_8_5x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p98_8_5y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p98_8_6(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p98_8_6x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p98_8_6y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p98_8_7(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p98_8_7x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p98_8_7y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p98_8_8(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p98_8_8x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p98_8_8y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p98_9_2(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p98_9_2x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p98_9_2y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p98_9_3(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p98_9_3x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p98_9_3y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p98_9_4(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p98_9_4x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p98_9_4y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p98_9_5(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p98_9_5x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p98_9_5y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p98_9_6(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p98_9_6x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p98_9_6y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p98_9_7(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p98_9_7x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p98_9_7y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p98_9_8(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p98_9_8x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p98_9_8y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p8_8(y);
}

/* order: 9, 9 */
static double eigen_quad_p99_2_2(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p99_2_2x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p99_2_2y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p99_2_3(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p99_2_3x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p99_2_3y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p99_2_4(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p99_2_4x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p99_2_4y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p99_2_5(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p99_2_5x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p99_2_5y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p99_2_6(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p99_2_6x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p99_2_6y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p99_2_7(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p99_2_7x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p99_2_7y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p99_2_8(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p99_2_8x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p99_2_8y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p99_2_9(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p99_2_9x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p99_2_9y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p99_3_2(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p99_3_2x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p99_3_2y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p99_3_3(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p99_3_3x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p99_3_3y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p99_3_4(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p99_3_4x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p99_3_4y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p99_3_5(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p99_3_5x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p99_3_5y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p99_3_6(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p99_3_6x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p99_3_6y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p99_3_7(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p99_3_7x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p99_3_7y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p99_3_8(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p99_3_8x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p99_3_8y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p99_3_9(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p99_3_9x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p99_3_9y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p99_4_2(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p99_4_2x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p99_4_2y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p99_4_3(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p99_4_3x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p99_4_3y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p99_4_4(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p99_4_4x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p99_4_4y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p99_4_5(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p99_4_5x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p99_4_5y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p99_4_6(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p99_4_6x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p99_4_6y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p99_4_7(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p99_4_7x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p99_4_7y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p99_4_8(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p99_4_8x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p99_4_8y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p99_4_9(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p99_4_9x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p99_4_9y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p99_5_2(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p99_5_2x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p99_5_2y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p99_5_3(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p99_5_3x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p99_5_3y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p99_5_4(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p99_5_4x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p99_5_4y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p99_5_5(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p99_5_5x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p99_5_5y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p99_5_6(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p99_5_6x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p99_5_6y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p99_5_7(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p99_5_7x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p99_5_7y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p99_5_8(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p99_5_8x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p99_5_8y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p99_5_9(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p99_5_9x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p99_5_9y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p99_6_2(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p99_6_2x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p99_6_2y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p99_6_3(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p99_6_3x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p99_6_3y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p99_6_4(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p99_6_4x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p99_6_4y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p99_6_5(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p99_6_5x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p99_6_5y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p99_6_6(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p99_6_6x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p99_6_6y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p99_6_7(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p99_6_7x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p99_6_7y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p99_6_8(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p99_6_8x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p99_6_8y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p99_6_9(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p99_6_9x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p99_6_9y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p99_7_2(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p99_7_2x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p99_7_2y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p99_7_3(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p99_7_3x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p99_7_3y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p99_7_4(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p99_7_4x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p99_7_4y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p99_7_5(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p99_7_5x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p99_7_5y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p99_7_6(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p99_7_6x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p99_7_6y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p99_7_7(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p99_7_7x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p99_7_7y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p99_7_8(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p99_7_8x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p99_7_8y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p99_7_9(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p99_7_9x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p99_7_9y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p99_8_2(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p99_8_2x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p99_8_2y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p99_8_3(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p99_8_3x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p99_8_3y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p99_8_4(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p99_8_4x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p99_8_4y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p99_8_5(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p99_8_5x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p99_8_5y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p99_8_6(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p99_8_6x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p99_8_6y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p99_8_7(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p99_8_7x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p99_8_7y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p99_8_8(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p99_8_8x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p99_8_8y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p99_8_9(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p99_8_9x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p99_8_9y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p99_9_2(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p99_9_2x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p99_9_2y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p99_9_3(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p99_9_3x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p99_9_3y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p99_9_4(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p99_9_4x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p99_9_4y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p99_9_5(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p99_9_5x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p99_9_5y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p99_9_6(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p99_9_6x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p99_9_6y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p99_9_7(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p99_9_7x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p99_9_7y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p99_9_8(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p99_9_8x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p99_9_8y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p99_9_9(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p99_9_9x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p99_9_9y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p9_9(y);
}

/* order: 9, 10 */
static double eigen_quad_p910_2_2(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p910_2_2x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p910_2_2y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p910_2_3(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p910_2_3x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p910_2_3y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p910_2_4(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p910_2_4x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p910_2_4y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p910_2_5(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p910_2_5x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p910_2_5y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p910_2_6(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p910_2_6x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p910_2_6y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p910_2_7(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p910_2_7x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p910_2_7y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p910_2_8(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p910_2_8x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p910_2_8y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p910_2_9(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p910_2_9x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p910_2_9y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p910_2_10(double x, double y)
{
  return  eigen_laplace_p9_2(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p910_2_10x(double x, double y)
{
  return  deigen_laplace_p9_2(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p910_2_10y(double x, double y)
{
  return  eigen_laplace_p9_2(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p910_3_2(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p910_3_2x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p910_3_2y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p910_3_3(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p910_3_3x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p910_3_3y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p910_3_4(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p910_3_4x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p910_3_4y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p910_3_5(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p910_3_5x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p910_3_5y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p910_3_6(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p910_3_6x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p910_3_6y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p910_3_7(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p910_3_7x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p910_3_7y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p910_3_8(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p910_3_8x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p910_3_8y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p910_3_9(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p910_3_9x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p910_3_9y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p910_3_10(double x, double y)
{
  return  eigen_laplace_p9_3(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p910_3_10x(double x, double y)
{
  return  deigen_laplace_p9_3(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p910_3_10y(double x, double y)
{
  return  eigen_laplace_p9_3(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p910_4_2(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p910_4_2x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p910_4_2y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p910_4_3(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p910_4_3x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p910_4_3y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p910_4_4(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p910_4_4x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p910_4_4y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p910_4_5(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p910_4_5x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p910_4_5y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p910_4_6(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p910_4_6x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p910_4_6y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p910_4_7(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p910_4_7x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p910_4_7y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p910_4_8(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p910_4_8x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p910_4_8y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p910_4_9(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p910_4_9x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p910_4_9y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p910_4_10(double x, double y)
{
  return  eigen_laplace_p9_4(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p910_4_10x(double x, double y)
{
  return  deigen_laplace_p9_4(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p910_4_10y(double x, double y)
{
  return  eigen_laplace_p9_4(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p910_5_2(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p910_5_2x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p910_5_2y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p910_5_3(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p910_5_3x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p910_5_3y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p910_5_4(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p910_5_4x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p910_5_4y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p910_5_5(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p910_5_5x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p910_5_5y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p910_5_6(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p910_5_6x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p910_5_6y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p910_5_7(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p910_5_7x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p910_5_7y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p910_5_8(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p910_5_8x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p910_5_8y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p910_5_9(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p910_5_9x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p910_5_9y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p910_5_10(double x, double y)
{
  return  eigen_laplace_p9_5(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p910_5_10x(double x, double y)
{
  return  deigen_laplace_p9_5(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p910_5_10y(double x, double y)
{
  return  eigen_laplace_p9_5(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p910_6_2(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p910_6_2x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p910_6_2y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p910_6_3(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p910_6_3x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p910_6_3y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p910_6_4(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p910_6_4x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p910_6_4y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p910_6_5(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p910_6_5x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p910_6_5y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p910_6_6(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p910_6_6x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p910_6_6y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p910_6_7(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p910_6_7x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p910_6_7y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p910_6_8(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p910_6_8x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p910_6_8y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p910_6_9(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p910_6_9x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p910_6_9y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p910_6_10(double x, double y)
{
  return  eigen_laplace_p9_6(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p910_6_10x(double x, double y)
{
  return  deigen_laplace_p9_6(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p910_6_10y(double x, double y)
{
  return  eigen_laplace_p9_6(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p910_7_2(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p910_7_2x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p910_7_2y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p910_7_3(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p910_7_3x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p910_7_3y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p910_7_4(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p910_7_4x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p910_7_4y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p910_7_5(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p910_7_5x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p910_7_5y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p910_7_6(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p910_7_6x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p910_7_6y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p910_7_7(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p910_7_7x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p910_7_7y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p910_7_8(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p910_7_8x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p910_7_8y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p910_7_9(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p910_7_9x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p910_7_9y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p910_7_10(double x, double y)
{
  return  eigen_laplace_p9_7(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p910_7_10x(double x, double y)
{
  return  deigen_laplace_p9_7(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p910_7_10y(double x, double y)
{
  return  eigen_laplace_p9_7(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p910_8_2(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p910_8_2x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p910_8_2y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p910_8_3(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p910_8_3x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p910_8_3y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p910_8_4(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p910_8_4x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p910_8_4y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p910_8_5(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p910_8_5x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p910_8_5y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p910_8_6(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p910_8_6x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p910_8_6y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p910_8_7(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p910_8_7x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p910_8_7y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p910_8_8(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p910_8_8x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p910_8_8y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p910_8_9(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p910_8_9x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p910_8_9y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p910_8_10(double x, double y)
{
  return  eigen_laplace_p9_8(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p910_8_10x(double x, double y)
{
  return  deigen_laplace_p9_8(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p910_8_10y(double x, double y)
{
  return  eigen_laplace_p9_8(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p910_9_2(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p910_9_2x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p910_9_2y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p910_9_3(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p910_9_3x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p910_9_3y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p910_9_4(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p910_9_4x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p910_9_4y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p910_9_5(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p910_9_5x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p910_9_5y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p910_9_6(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p910_9_6x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p910_9_6y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p910_9_7(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p910_9_7x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p910_9_7y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p910_9_8(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p910_9_8x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p910_9_8y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p910_9_9(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p910_9_9x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p910_9_9y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p910_9_10(double x, double y)
{
  return  eigen_laplace_p9_9(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p910_9_10x(double x, double y)
{
  return  deigen_laplace_p9_9(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p910_9_10y(double x, double y)
{
  return  eigen_laplace_p9_9(x) * deigen_laplace_p10_10(y);
}

/* order: 10, 2 */
static double eigen_quad_p102_2_2(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p102_2_2x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p102_2_2y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p102_3_2(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p102_3_2x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p102_3_2y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p102_4_2(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p102_4_2x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p102_4_2y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p102_5_2(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p102_5_2x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p102_5_2y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p102_6_2(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p102_6_2x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p102_6_2y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p102_7_2(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p102_7_2x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p102_7_2y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p102_8_2(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p102_8_2x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p102_8_2y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p102_9_2(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p102_9_2x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p102_9_2y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p2_2(y);
}

static double eigen_quad_p102_10_2(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p102_10_2x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p2_2(y);
}

static double eigen_quad_p102_10_2y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p2_2(y);
}

/* order: 10, 3 */
static double eigen_quad_p103_2_2(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p103_2_2x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p103_2_2y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p103_2_3(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p103_2_3x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p103_2_3y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p103_3_2(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p103_3_2x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p103_3_2y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p103_3_3(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p103_3_3x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p103_3_3y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p103_4_2(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p103_4_2x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p103_4_2y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p103_4_3(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p103_4_3x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p103_4_3y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p103_5_2(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p103_5_2x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p103_5_2y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p103_5_3(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p103_5_3x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p103_5_3y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p103_6_2(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p103_6_2x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p103_6_2y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p103_6_3(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p103_6_3x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p103_6_3y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p103_7_2(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p103_7_2x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p103_7_2y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p103_7_3(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p103_7_3x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p103_7_3y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p103_8_2(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p103_8_2x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p103_8_2y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p103_8_3(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p103_8_3x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p103_8_3y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p103_9_2(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p103_9_2x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p103_9_2y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p103_9_3(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p103_9_3x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p103_9_3y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p3_3(y);
}

static double eigen_quad_p103_10_2(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p103_10_2x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p3_2(y);
}

static double eigen_quad_p103_10_2y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p3_2(y);
}

static double eigen_quad_p103_10_3(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p103_10_3x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p3_3(y);
}

static double eigen_quad_p103_10_3y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p3_3(y);
}

/* order: 10, 4 */
static double eigen_quad_p104_2_2(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p104_2_2x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p104_2_2y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p104_2_3(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p104_2_3x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p104_2_3y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p104_2_4(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p104_2_4x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p104_2_4y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p104_3_2(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p104_3_2x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p104_3_2y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p104_3_3(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p104_3_3x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p104_3_3y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p104_3_4(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p104_3_4x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p104_3_4y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p104_4_2(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p104_4_2x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p104_4_2y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p104_4_3(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p104_4_3x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p104_4_3y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p104_4_4(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p104_4_4x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p104_4_4y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p104_5_2(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p104_5_2x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p104_5_2y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p104_5_3(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p104_5_3x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p104_5_3y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p104_5_4(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p104_5_4x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p104_5_4y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p104_6_2(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p104_6_2x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p104_6_2y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p104_6_3(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p104_6_3x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p104_6_3y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p104_6_4(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p104_6_4x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p104_6_4y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p104_7_2(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p104_7_2x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p104_7_2y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p104_7_3(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p104_7_3x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p104_7_3y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p104_7_4(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p104_7_4x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p104_7_4y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p104_8_2(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p104_8_2x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p104_8_2y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p104_8_3(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p104_8_3x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p104_8_3y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p104_8_4(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p104_8_4x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p104_8_4y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p104_9_2(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p104_9_2x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p104_9_2y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p104_9_3(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p104_9_3x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p104_9_3y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p104_9_4(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p104_9_4x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p104_9_4y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p4_4(y);
}

static double eigen_quad_p104_10_2(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p104_10_2x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p4_2(y);
}

static double eigen_quad_p104_10_2y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p4_2(y);
}

static double eigen_quad_p104_10_3(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p104_10_3x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p4_3(y);
}

static double eigen_quad_p104_10_3y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p4_3(y);
}

static double eigen_quad_p104_10_4(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p104_10_4x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p4_4(y);
}

static double eigen_quad_p104_10_4y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p4_4(y);
}

/* order: 10, 5 */
static double eigen_quad_p105_2_2(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p105_2_2x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p105_2_2y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p105_2_3(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p105_2_3x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p105_2_3y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p105_2_4(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p105_2_4x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p105_2_4y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p105_2_5(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p105_2_5x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p105_2_5y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p105_3_2(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p105_3_2x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p105_3_2y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p105_3_3(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p105_3_3x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p105_3_3y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p105_3_4(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p105_3_4x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p105_3_4y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p105_3_5(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p105_3_5x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p105_3_5y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p105_4_2(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p105_4_2x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p105_4_2y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p105_4_3(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p105_4_3x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p105_4_3y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p105_4_4(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p105_4_4x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p105_4_4y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p105_4_5(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p105_4_5x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p105_4_5y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p105_5_2(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p105_5_2x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p105_5_2y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p105_5_3(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p105_5_3x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p105_5_3y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p105_5_4(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p105_5_4x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p105_5_4y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p105_5_5(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p105_5_5x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p105_5_5y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p105_6_2(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p105_6_2x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p105_6_2y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p105_6_3(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p105_6_3x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p105_6_3y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p105_6_4(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p105_6_4x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p105_6_4y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p105_6_5(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p105_6_5x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p105_6_5y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p105_7_2(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p105_7_2x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p105_7_2y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p105_7_3(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p105_7_3x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p105_7_3y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p105_7_4(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p105_7_4x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p105_7_4y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p105_7_5(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p105_7_5x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p105_7_5y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p105_8_2(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p105_8_2x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p105_8_2y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p105_8_3(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p105_8_3x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p105_8_3y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p105_8_4(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p105_8_4x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p105_8_4y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p105_8_5(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p105_8_5x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p105_8_5y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p105_9_2(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p105_9_2x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p105_9_2y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p105_9_3(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p105_9_3x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p105_9_3y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p105_9_4(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p105_9_4x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p105_9_4y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p105_9_5(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p105_9_5x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p105_9_5y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p5_5(y);
}

static double eigen_quad_p105_10_2(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p105_10_2x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p5_2(y);
}

static double eigen_quad_p105_10_2y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p5_2(y);
}

static double eigen_quad_p105_10_3(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p105_10_3x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p5_3(y);
}

static double eigen_quad_p105_10_3y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p5_3(y);
}

static double eigen_quad_p105_10_4(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p105_10_4x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p5_4(y);
}

static double eigen_quad_p105_10_4y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p5_4(y);
}

static double eigen_quad_p105_10_5(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p105_10_5x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p5_5(y);
}

static double eigen_quad_p105_10_5y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p5_5(y);
}

/* order: 10, 6 */
static double eigen_quad_p106_2_2(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p106_2_2x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p106_2_2y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p106_2_3(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p106_2_3x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p106_2_3y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p106_2_4(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p106_2_4x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p106_2_4y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p106_2_5(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p106_2_5x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p106_2_5y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p106_2_6(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p106_2_6x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p106_2_6y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p106_3_2(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p106_3_2x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p106_3_2y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p106_3_3(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p106_3_3x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p106_3_3y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p106_3_4(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p106_3_4x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p106_3_4y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p106_3_5(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p106_3_5x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p106_3_5y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p106_3_6(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p106_3_6x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p106_3_6y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p106_4_2(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p106_4_2x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p106_4_2y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p106_4_3(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p106_4_3x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p106_4_3y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p106_4_4(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p106_4_4x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p106_4_4y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p106_4_5(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p106_4_5x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p106_4_5y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p106_4_6(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p106_4_6x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p106_4_6y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p106_5_2(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p106_5_2x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p106_5_2y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p106_5_3(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p106_5_3x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p106_5_3y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p106_5_4(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p106_5_4x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p106_5_4y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p106_5_5(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p106_5_5x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p106_5_5y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p106_5_6(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p106_5_6x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p106_5_6y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p106_6_2(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p106_6_2x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p106_6_2y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p106_6_3(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p106_6_3x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p106_6_3y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p106_6_4(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p106_6_4x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p106_6_4y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p106_6_5(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p106_6_5x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p106_6_5y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p106_6_6(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p106_6_6x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p106_6_6y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p106_7_2(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p106_7_2x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p106_7_2y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p106_7_3(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p106_7_3x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p106_7_3y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p106_7_4(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p106_7_4x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p106_7_4y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p106_7_5(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p106_7_5x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p106_7_5y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p106_7_6(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p106_7_6x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p106_7_6y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p106_8_2(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p106_8_2x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p106_8_2y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p106_8_3(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p106_8_3x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p106_8_3y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p106_8_4(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p106_8_4x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p106_8_4y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p106_8_5(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p106_8_5x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p106_8_5y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p106_8_6(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p106_8_6x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p106_8_6y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p106_9_2(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p106_9_2x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p106_9_2y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p106_9_3(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p106_9_3x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p106_9_3y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p106_9_4(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p106_9_4x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p106_9_4y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p106_9_5(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p106_9_5x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p106_9_5y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p106_9_6(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p106_9_6x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p106_9_6y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p6_6(y);
}

static double eigen_quad_p106_10_2(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p106_10_2x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p6_2(y);
}

static double eigen_quad_p106_10_2y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p6_2(y);
}

static double eigen_quad_p106_10_3(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p106_10_3x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p6_3(y);
}

static double eigen_quad_p106_10_3y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p6_3(y);
}

static double eigen_quad_p106_10_4(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p106_10_4x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p6_4(y);
}

static double eigen_quad_p106_10_4y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p6_4(y);
}

static double eigen_quad_p106_10_5(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p106_10_5x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p6_5(y);
}

static double eigen_quad_p106_10_5y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p6_5(y);
}

static double eigen_quad_p106_10_6(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p106_10_6x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p6_6(y);
}

static double eigen_quad_p106_10_6y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p6_6(y);
}

/* order: 10, 7 */
static double eigen_quad_p107_2_2(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p107_2_2x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p107_2_2y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p107_2_3(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p107_2_3x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p107_2_3y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p107_2_4(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p107_2_4x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p107_2_4y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p107_2_5(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p107_2_5x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p107_2_5y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p107_2_6(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p107_2_6x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p107_2_6y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p107_2_7(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p107_2_7x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p107_2_7y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p107_3_2(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p107_3_2x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p107_3_2y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p107_3_3(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p107_3_3x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p107_3_3y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p107_3_4(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p107_3_4x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p107_3_4y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p107_3_5(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p107_3_5x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p107_3_5y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p107_3_6(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p107_3_6x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p107_3_6y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p107_3_7(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p107_3_7x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p107_3_7y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p107_4_2(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p107_4_2x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p107_4_2y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p107_4_3(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p107_4_3x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p107_4_3y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p107_4_4(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p107_4_4x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p107_4_4y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p107_4_5(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p107_4_5x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p107_4_5y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p107_4_6(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p107_4_6x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p107_4_6y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p107_4_7(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p107_4_7x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p107_4_7y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p107_5_2(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p107_5_2x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p107_5_2y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p107_5_3(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p107_5_3x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p107_5_3y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p107_5_4(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p107_5_4x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p107_5_4y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p107_5_5(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p107_5_5x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p107_5_5y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p107_5_6(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p107_5_6x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p107_5_6y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p107_5_7(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p107_5_7x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p107_5_7y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p107_6_2(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p107_6_2x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p107_6_2y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p107_6_3(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p107_6_3x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p107_6_3y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p107_6_4(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p107_6_4x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p107_6_4y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p107_6_5(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p107_6_5x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p107_6_5y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p107_6_6(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p107_6_6x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p107_6_6y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p107_6_7(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p107_6_7x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p107_6_7y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p107_7_2(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p107_7_2x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p107_7_2y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p107_7_3(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p107_7_3x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p107_7_3y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p107_7_4(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p107_7_4x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p107_7_4y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p107_7_5(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p107_7_5x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p107_7_5y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p107_7_6(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p107_7_6x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p107_7_6y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p107_7_7(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p107_7_7x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p107_7_7y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p107_8_2(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p107_8_2x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p107_8_2y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p107_8_3(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p107_8_3x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p107_8_3y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p107_8_4(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p107_8_4x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p107_8_4y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p107_8_5(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p107_8_5x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p107_8_5y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p107_8_6(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p107_8_6x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p107_8_6y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p107_8_7(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p107_8_7x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p107_8_7y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p107_9_2(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p107_9_2x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p107_9_2y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p107_9_3(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p107_9_3x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p107_9_3y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p107_9_4(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p107_9_4x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p107_9_4y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p107_9_5(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p107_9_5x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p107_9_5y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p107_9_6(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p107_9_6x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p107_9_6y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p107_9_7(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p107_9_7x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p107_9_7y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p7_7(y);
}

static double eigen_quad_p107_10_2(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p107_10_2x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p7_2(y);
}

static double eigen_quad_p107_10_2y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p7_2(y);
}

static double eigen_quad_p107_10_3(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p107_10_3x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p7_3(y);
}

static double eigen_quad_p107_10_3y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p7_3(y);
}

static double eigen_quad_p107_10_4(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p107_10_4x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p7_4(y);
}

static double eigen_quad_p107_10_4y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p7_4(y);
}

static double eigen_quad_p107_10_5(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p107_10_5x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p7_5(y);
}

static double eigen_quad_p107_10_5y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p7_5(y);
}

static double eigen_quad_p107_10_6(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p107_10_6x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p7_6(y);
}

static double eigen_quad_p107_10_6y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p7_6(y);
}

static double eigen_quad_p107_10_7(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p107_10_7x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p7_7(y);
}

static double eigen_quad_p107_10_7y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p7_7(y);
}

/* order: 10, 8 */
static double eigen_quad_p108_2_2(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p108_2_2x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p108_2_2y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p108_2_3(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p108_2_3x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p108_2_3y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p108_2_4(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p108_2_4x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p108_2_4y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p108_2_5(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p108_2_5x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p108_2_5y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p108_2_6(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p108_2_6x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p108_2_6y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p108_2_7(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p108_2_7x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p108_2_7y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p108_2_8(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p108_2_8x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p108_2_8y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p108_3_2(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p108_3_2x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p108_3_2y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p108_3_3(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p108_3_3x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p108_3_3y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p108_3_4(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p108_3_4x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p108_3_4y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p108_3_5(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p108_3_5x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p108_3_5y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p108_3_6(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p108_3_6x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p108_3_6y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p108_3_7(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p108_3_7x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p108_3_7y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p108_3_8(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p108_3_8x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p108_3_8y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p108_4_2(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p108_4_2x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p108_4_2y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p108_4_3(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p108_4_3x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p108_4_3y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p108_4_4(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p108_4_4x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p108_4_4y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p108_4_5(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p108_4_5x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p108_4_5y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p108_4_6(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p108_4_6x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p108_4_6y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p108_4_7(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p108_4_7x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p108_4_7y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p108_4_8(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p108_4_8x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p108_4_8y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p108_5_2(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p108_5_2x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p108_5_2y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p108_5_3(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p108_5_3x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p108_5_3y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p108_5_4(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p108_5_4x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p108_5_4y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p108_5_5(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p108_5_5x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p108_5_5y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p108_5_6(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p108_5_6x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p108_5_6y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p108_5_7(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p108_5_7x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p108_5_7y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p108_5_8(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p108_5_8x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p108_5_8y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p108_6_2(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p108_6_2x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p108_6_2y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p108_6_3(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p108_6_3x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p108_6_3y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p108_6_4(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p108_6_4x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p108_6_4y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p108_6_5(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p108_6_5x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p108_6_5y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p108_6_6(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p108_6_6x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p108_6_6y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p108_6_7(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p108_6_7x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p108_6_7y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p108_6_8(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p108_6_8x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p108_6_8y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p108_7_2(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p108_7_2x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p108_7_2y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p108_7_3(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p108_7_3x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p108_7_3y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p108_7_4(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p108_7_4x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p108_7_4y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p108_7_5(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p108_7_5x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p108_7_5y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p108_7_6(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p108_7_6x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p108_7_6y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p108_7_7(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p108_7_7x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p108_7_7y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p108_7_8(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p108_7_8x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p108_7_8y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p108_8_2(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p108_8_2x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p108_8_2y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p108_8_3(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p108_8_3x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p108_8_3y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p108_8_4(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p108_8_4x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p108_8_4y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p108_8_5(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p108_8_5x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p108_8_5y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p108_8_6(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p108_8_6x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p108_8_6y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p108_8_7(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p108_8_7x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p108_8_7y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p108_8_8(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p108_8_8x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p108_8_8y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p108_9_2(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p108_9_2x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p108_9_2y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p108_9_3(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p108_9_3x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p108_9_3y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p108_9_4(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p108_9_4x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p108_9_4y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p108_9_5(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p108_9_5x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p108_9_5y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p108_9_6(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p108_9_6x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p108_9_6y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p108_9_7(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p108_9_7x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p108_9_7y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p108_9_8(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p108_9_8x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p108_9_8y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p8_8(y);
}

static double eigen_quad_p108_10_2(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p108_10_2x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p8_2(y);
}

static double eigen_quad_p108_10_2y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p8_2(y);
}

static double eigen_quad_p108_10_3(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p108_10_3x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p8_3(y);
}

static double eigen_quad_p108_10_3y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p8_3(y);
}

static double eigen_quad_p108_10_4(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p108_10_4x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p8_4(y);
}

static double eigen_quad_p108_10_4y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p8_4(y);
}

static double eigen_quad_p108_10_5(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p108_10_5x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p8_5(y);
}

static double eigen_quad_p108_10_5y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p8_5(y);
}

static double eigen_quad_p108_10_6(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p108_10_6x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p8_6(y);
}

static double eigen_quad_p108_10_6y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p8_6(y);
}

static double eigen_quad_p108_10_7(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p108_10_7x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p8_7(y);
}

static double eigen_quad_p108_10_7y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p8_7(y);
}

static double eigen_quad_p108_10_8(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p108_10_8x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p8_8(y);
}

static double eigen_quad_p108_10_8y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p8_8(y);
}

/* order: 10, 9 */
static double eigen_quad_p109_2_2(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p109_2_2x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p109_2_2y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p109_2_3(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p109_2_3x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p109_2_3y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p109_2_4(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p109_2_4x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p109_2_4y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p109_2_5(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p109_2_5x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p109_2_5y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p109_2_6(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p109_2_6x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p109_2_6y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p109_2_7(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p109_2_7x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p109_2_7y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p109_2_8(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p109_2_8x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p109_2_8y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p109_2_9(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p109_2_9x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p109_2_9y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p109_3_2(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p109_3_2x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p109_3_2y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p109_3_3(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p109_3_3x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p109_3_3y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p109_3_4(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p109_3_4x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p109_3_4y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p109_3_5(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p109_3_5x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p109_3_5y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p109_3_6(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p109_3_6x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p109_3_6y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p109_3_7(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p109_3_7x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p109_3_7y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p109_3_8(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p109_3_8x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p109_3_8y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p109_3_9(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p109_3_9x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p109_3_9y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p109_4_2(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p109_4_2x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p109_4_2y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p109_4_3(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p109_4_3x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p109_4_3y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p109_4_4(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p109_4_4x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p109_4_4y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p109_4_5(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p109_4_5x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p109_4_5y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p109_4_6(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p109_4_6x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p109_4_6y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p109_4_7(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p109_4_7x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p109_4_7y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p109_4_8(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p109_4_8x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p109_4_8y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p109_4_9(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p109_4_9x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p109_4_9y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p109_5_2(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p109_5_2x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p109_5_2y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p109_5_3(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p109_5_3x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p109_5_3y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p109_5_4(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p109_5_4x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p109_5_4y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p109_5_5(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p109_5_5x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p109_5_5y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p109_5_6(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p109_5_6x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p109_5_6y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p109_5_7(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p109_5_7x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p109_5_7y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p109_5_8(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p109_5_8x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p109_5_8y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p109_5_9(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p109_5_9x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p109_5_9y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p109_6_2(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p109_6_2x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p109_6_2y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p109_6_3(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p109_6_3x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p109_6_3y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p109_6_4(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p109_6_4x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p109_6_4y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p109_6_5(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p109_6_5x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p109_6_5y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p109_6_6(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p109_6_6x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p109_6_6y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p109_6_7(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p109_6_7x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p109_6_7y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p109_6_8(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p109_6_8x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p109_6_8y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p109_6_9(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p109_6_9x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p109_6_9y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p109_7_2(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p109_7_2x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p109_7_2y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p109_7_3(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p109_7_3x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p109_7_3y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p109_7_4(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p109_7_4x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p109_7_4y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p109_7_5(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p109_7_5x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p109_7_5y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p109_7_6(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p109_7_6x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p109_7_6y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p109_7_7(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p109_7_7x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p109_7_7y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p109_7_8(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p109_7_8x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p109_7_8y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p109_7_9(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p109_7_9x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p109_7_9y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p109_8_2(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p109_8_2x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p109_8_2y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p109_8_3(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p109_8_3x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p109_8_3y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p109_8_4(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p109_8_4x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p109_8_4y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p109_8_5(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p109_8_5x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p109_8_5y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p109_8_6(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p109_8_6x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p109_8_6y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p109_8_7(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p109_8_7x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p109_8_7y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p109_8_8(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p109_8_8x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p109_8_8y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p109_8_9(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p109_8_9x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p109_8_9y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p109_9_2(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p109_9_2x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p109_9_2y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p109_9_3(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p109_9_3x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p109_9_3y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p109_9_4(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p109_9_4x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p109_9_4y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p109_9_5(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p109_9_5x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p109_9_5y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p109_9_6(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p109_9_6x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p109_9_6y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p109_9_7(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p109_9_7x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p109_9_7y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p109_9_8(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p109_9_8x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p109_9_8y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p109_9_9(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p109_9_9x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p109_9_9y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p9_9(y);
}

static double eigen_quad_p109_10_2(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p109_10_2x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p9_2(y);
}

static double eigen_quad_p109_10_2y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p9_2(y);
}

static double eigen_quad_p109_10_3(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p109_10_3x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p9_3(y);
}

static double eigen_quad_p109_10_3y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p9_3(y);
}

static double eigen_quad_p109_10_4(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p109_10_4x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p9_4(y);
}

static double eigen_quad_p109_10_4y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p9_4(y);
}

static double eigen_quad_p109_10_5(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p109_10_5x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p9_5(y);
}

static double eigen_quad_p109_10_5y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p9_5(y);
}

static double eigen_quad_p109_10_6(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p109_10_6x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p9_6(y);
}

static double eigen_quad_p109_10_6y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p9_6(y);
}

static double eigen_quad_p109_10_7(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p109_10_7x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p9_7(y);
}

static double eigen_quad_p109_10_7y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p9_7(y);
}

static double eigen_quad_p109_10_8(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p109_10_8x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p9_8(y);
}

static double eigen_quad_p109_10_8y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p9_8(y);
}

static double eigen_quad_p109_10_9(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p109_10_9x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p9_9(y);
}

static double eigen_quad_p109_10_9y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p9_9(y);
}

/* order: 10, 10 */
static double eigen_quad_p1010_2_2(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p1010_2_2x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p1010_2_2y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p1010_2_3(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p1010_2_3x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p1010_2_3y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p1010_2_4(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p1010_2_4x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p1010_2_4y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p1010_2_5(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p1010_2_5x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p1010_2_5y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p1010_2_6(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p1010_2_6x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p1010_2_6y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p1010_2_7(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p1010_2_7x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p1010_2_7y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p1010_2_8(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p1010_2_8x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p1010_2_8y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p1010_2_9(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p1010_2_9x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p1010_2_9y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p1010_2_10(double x, double y)
{
  return  eigen_laplace_p10_2(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p1010_2_10x(double x, double y)
{
  return  deigen_laplace_p10_2(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p1010_2_10y(double x, double y)
{
  return  eigen_laplace_p10_2(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p1010_3_2(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p1010_3_2x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p1010_3_2y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p1010_3_3(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p1010_3_3x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p1010_3_3y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p1010_3_4(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p1010_3_4x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p1010_3_4y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p1010_3_5(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p1010_3_5x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p1010_3_5y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p1010_3_6(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p1010_3_6x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p1010_3_6y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p1010_3_7(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p1010_3_7x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p1010_3_7y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p1010_3_8(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p1010_3_8x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p1010_3_8y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p1010_3_9(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p1010_3_9x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p1010_3_9y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p1010_3_10(double x, double y)
{
  return  eigen_laplace_p10_3(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p1010_3_10x(double x, double y)
{
  return  deigen_laplace_p10_3(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p1010_3_10y(double x, double y)
{
  return  eigen_laplace_p10_3(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p1010_4_2(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p1010_4_2x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p1010_4_2y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p1010_4_3(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p1010_4_3x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p1010_4_3y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p1010_4_4(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p1010_4_4x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p1010_4_4y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p1010_4_5(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p1010_4_5x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p1010_4_5y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p1010_4_6(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p1010_4_6x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p1010_4_6y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p1010_4_7(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p1010_4_7x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p1010_4_7y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p1010_4_8(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p1010_4_8x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p1010_4_8y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p1010_4_9(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p1010_4_9x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p1010_4_9y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p1010_4_10(double x, double y)
{
  return  eigen_laplace_p10_4(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p1010_4_10x(double x, double y)
{
  return  deigen_laplace_p10_4(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p1010_4_10y(double x, double y)
{
  return  eigen_laplace_p10_4(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p1010_5_2(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p1010_5_2x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p1010_5_2y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p1010_5_3(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p1010_5_3x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p1010_5_3y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p1010_5_4(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p1010_5_4x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p1010_5_4y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p1010_5_5(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p1010_5_5x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p1010_5_5y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p1010_5_6(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p1010_5_6x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p1010_5_6y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p1010_5_7(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p1010_5_7x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p1010_5_7y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p1010_5_8(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p1010_5_8x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p1010_5_8y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p1010_5_9(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p1010_5_9x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p1010_5_9y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p1010_5_10(double x, double y)
{
  return  eigen_laplace_p10_5(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p1010_5_10x(double x, double y)
{
  return  deigen_laplace_p10_5(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p1010_5_10y(double x, double y)
{
  return  eigen_laplace_p10_5(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p1010_6_2(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p1010_6_2x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p1010_6_2y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p1010_6_3(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p1010_6_3x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p1010_6_3y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p1010_6_4(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p1010_6_4x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p1010_6_4y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p1010_6_5(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p1010_6_5x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p1010_6_5y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p1010_6_6(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p1010_6_6x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p1010_6_6y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p1010_6_7(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p1010_6_7x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p1010_6_7y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p1010_6_8(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p1010_6_8x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p1010_6_8y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p1010_6_9(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p1010_6_9x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p1010_6_9y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p1010_6_10(double x, double y)
{
  return  eigen_laplace_p10_6(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p1010_6_10x(double x, double y)
{
  return  deigen_laplace_p10_6(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p1010_6_10y(double x, double y)
{
  return  eigen_laplace_p10_6(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p1010_7_2(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p1010_7_2x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p1010_7_2y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p1010_7_3(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p1010_7_3x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p1010_7_3y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p1010_7_4(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p1010_7_4x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p1010_7_4y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p1010_7_5(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p1010_7_5x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p1010_7_5y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p1010_7_6(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p1010_7_6x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p1010_7_6y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p1010_7_7(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p1010_7_7x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p1010_7_7y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p1010_7_8(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p1010_7_8x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p1010_7_8y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p1010_7_9(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p1010_7_9x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p1010_7_9y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p1010_7_10(double x, double y)
{
  return  eigen_laplace_p10_7(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p1010_7_10x(double x, double y)
{
  return  deigen_laplace_p10_7(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p1010_7_10y(double x, double y)
{
  return  eigen_laplace_p10_7(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p1010_8_2(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p1010_8_2x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p1010_8_2y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p1010_8_3(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p1010_8_3x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p1010_8_3y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p1010_8_4(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p1010_8_4x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p1010_8_4y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p1010_8_5(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p1010_8_5x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p1010_8_5y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p1010_8_6(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p1010_8_6x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p1010_8_6y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p1010_8_7(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p1010_8_7x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p1010_8_7y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p1010_8_8(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p1010_8_8x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p1010_8_8y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p1010_8_9(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p1010_8_9x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p1010_8_9y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p1010_8_10(double x, double y)
{
  return  eigen_laplace_p10_8(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p1010_8_10x(double x, double y)
{
  return  deigen_laplace_p10_8(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p1010_8_10y(double x, double y)
{
  return  eigen_laplace_p10_8(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p1010_9_2(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p1010_9_2x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p1010_9_2y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p1010_9_3(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p1010_9_3x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p1010_9_3y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p1010_9_4(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p1010_9_4x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p1010_9_4y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p1010_9_5(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p1010_9_5x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p1010_9_5y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p1010_9_6(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p1010_9_6x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p1010_9_6y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p1010_9_7(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p1010_9_7x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p1010_9_7y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p1010_9_8(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p1010_9_8x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p1010_9_8y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p1010_9_9(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p1010_9_9x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p1010_9_9y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p1010_9_10(double x, double y)
{
  return  eigen_laplace_p10_9(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p1010_9_10x(double x, double y)
{
  return  deigen_laplace_p10_9(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p1010_9_10y(double x, double y)
{
  return  eigen_laplace_p10_9(x) * deigen_laplace_p10_10(y);
}

static double eigen_quad_p1010_10_2(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p1010_10_2x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p10_2(y);
}

static double eigen_quad_p1010_10_2y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p10_2(y);
}

static double eigen_quad_p1010_10_3(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p1010_10_3x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p10_3(y);
}

static double eigen_quad_p1010_10_3y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p10_3(y);
}

static double eigen_quad_p1010_10_4(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p1010_10_4x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p10_4(y);
}

static double eigen_quad_p1010_10_4y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p10_4(y);
}

static double eigen_quad_p1010_10_5(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p1010_10_5x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p10_5(y);
}

static double eigen_quad_p1010_10_5y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p10_5(y);
}

static double eigen_quad_p1010_10_6(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p1010_10_6x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p10_6(y);
}

static double eigen_quad_p1010_10_6y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p10_6(y);
}

static double eigen_quad_p1010_10_7(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p1010_10_7x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p10_7(y);
}

static double eigen_quad_p1010_10_7y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p10_7(y);
}

static double eigen_quad_p1010_10_8(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p1010_10_8x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p10_8(y);
}

static double eigen_quad_p1010_10_8y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p10_8(y);
}

static double eigen_quad_p1010_10_9(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p1010_10_9x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p10_9(y);
}

static double eigen_quad_p1010_10_9y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p10_9(y);
}

static double eigen_quad_p1010_10_10(double x, double y)
{
  return  eigen_laplace_p10_10(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p1010_10_10x(double x, double y)
{
  return  deigen_laplace_p10_10(x) * eigen_laplace_p10_10(y);
}

static double eigen_quad_p1010_10_10y(double x, double y)
{
  return  eigen_laplace_p10_10(x) * deigen_laplace_p10_10(y);
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////

static Shapeset::shape_fn_t eigen_quad_fn[] =
{
  eigen_quad_0_0,   eigen_quad_0_1,   eigen_quad_0_2,   eigen_quad_0_3_0, eigen_quad_0_3_1,  eigen_quad_0_4,   eigen_quad_0_5_0,
  eigen_quad_0_5_1, eigen_quad_0_6,   eigen_quad_0_7_0, eigen_quad_0_7_1, eigen_quad_0_8,   eigen_quad_0_9_0,  eigen_quad_0_9_1,
  eigen_quad_0_10,  eigen_quad_1_0,  eigen_quad_1_1,   eigen_quad_1_2,   eigen_quad_1_3_0, eigen_quad_1_3_1, eigen_quad_1_4,
     eigen_quad_1_5_0, eigen_quad_1_5_1, eigen_quad_1_6,  eigen_quad_1_7_0, eigen_quad_1_7_1, eigen_quad_1_8,   eigen_quad_1_9_0,
     eigen_quad_1_9_1,  eigen_quad_1_10,  eigen_quad_2_0,   eigen_quad_2_1,   eigen_quad_3_0_0, eigen_quad_3_0_1, eigen_quad_3_1_0,
     eigen_quad_3_1_1, eigen_quad_4_0,   eigen_quad_4_1,   eigen_quad_5_0_0, eigen_quad_5_0_1, eigen_quad_5_1_0, eigen_quad_5_1_1,
     eigen_quad_6_0,   eigen_quad_6_1,   eigen_quad_7_0_0, eigen_quad_7_0_1, eigen_quad_7_1_0, eigen_quad_7_1_1, eigen_quad_8_0,
     eigen_quad_8_1,   eigen_quad_9_0_0, eigen_quad_9_0_1, eigen_quad_9_1_0, eigen_quad_9_1_1, eigen_quad_10_0,  eigen_quad_10_1,

  eigen_quad_p22_2_2, eigen_quad_p23_2_2, eigen_quad_p23_2_3, eigen_quad_p24_2_2, eigen_quad_p24_2_3, eigen_quad_p24_2_4,
  eigen_quad_p25_2_2, eigen_quad_p25_2_3, eigen_quad_p25_2_4, eigen_quad_p25_2_5, eigen_quad_p26_2_2, eigen_quad_p26_2_3,
  eigen_quad_p26_2_4, eigen_quad_p26_2_5, eigen_quad_p26_2_6, eigen_quad_p27_2_2, eigen_quad_p27_2_3, eigen_quad_p27_2_4,
  eigen_quad_p27_2_5, eigen_quad_p27_2_6, eigen_quad_p27_2_7, eigen_quad_p28_2_2, eigen_quad_p28_2_3, eigen_quad_p28_2_4,
  eigen_quad_p28_2_5, eigen_quad_p28_2_6, eigen_quad_p28_2_7, eigen_quad_p28_2_8, eigen_quad_p29_2_2, eigen_quad_p29_2_3,
  eigen_quad_p29_2_4, eigen_quad_p29_2_5, eigen_quad_p29_2_6, eigen_quad_p29_2_7, eigen_quad_p29_2_8, eigen_quad_p29_2_9,
  eigen_quad_p210_2_2, eigen_quad_p210_2_3, eigen_quad_p210_2_4, eigen_quad_p210_2_5, eigen_quad_p210_2_6, eigen_quad_p210_2_7,
  eigen_quad_p210_2_8, eigen_quad_p210_2_9, eigen_quad_p210_2_10,

  eigen_quad_p32_2_2, eigen_quad_p32_3_2, eigen_quad_p33_2_2, eigen_quad_p33_2_3, eigen_quad_p33_3_2, eigen_quad_p33_3_3,
  eigen_quad_p34_2_2, eigen_quad_p34_2_3, eigen_quad_p34_2_4, eigen_quad_p34_3_2, eigen_quad_p34_3_3, eigen_quad_p34_3_4,
  eigen_quad_p35_2_2, eigen_quad_p35_2_3, eigen_quad_p35_2_4, eigen_quad_p35_2_5, eigen_quad_p35_3_2, eigen_quad_p35_3_3,
  eigen_quad_p35_3_4, eigen_quad_p35_3_5, eigen_quad_p36_2_2, eigen_quad_p36_2_3, eigen_quad_p36_2_4, eigen_quad_p36_2_5,
  eigen_quad_p36_2_6, eigen_quad_p36_3_2, eigen_quad_p36_3_3, eigen_quad_p36_3_4, eigen_quad_p36_3_5, eigen_quad_p36_3_6,
  eigen_quad_p37_2_2, eigen_quad_p37_2_3, eigen_quad_p37_2_4, eigen_quad_p37_2_5, eigen_quad_p37_2_6, eigen_quad_p37_2_7,
  eigen_quad_p37_3_2, eigen_quad_p37_3_3, eigen_quad_p37_3_4, eigen_quad_p37_3_5, eigen_quad_p37_3_6, eigen_quad_p37_3_7,
  eigen_quad_p38_2_2, eigen_quad_p38_2_3, eigen_quad_p38_2_4, eigen_quad_p38_2_5, eigen_quad_p38_2_6, eigen_quad_p38_2_7,
  eigen_quad_p38_2_8, eigen_quad_p38_3_2, eigen_quad_p38_3_3, eigen_quad_p38_3_4, eigen_quad_p38_3_5, eigen_quad_p38_3_6,
  eigen_quad_p38_3_7, eigen_quad_p38_3_8, eigen_quad_p39_2_2, eigen_quad_p39_2_3, eigen_quad_p39_2_4, eigen_quad_p39_2_5,
  eigen_quad_p39_2_6, eigen_quad_p39_2_7, eigen_quad_p39_2_8, eigen_quad_p39_2_9, eigen_quad_p39_3_2, eigen_quad_p39_3_3,
  eigen_quad_p39_3_4, eigen_quad_p39_3_5, eigen_quad_p39_3_6, eigen_quad_p39_3_7, eigen_quad_p39_3_8, eigen_quad_p39_3_9,
  eigen_quad_p310_2_2, eigen_quad_p310_2_3, eigen_quad_p310_2_4, eigen_quad_p310_2_5, eigen_quad_p310_2_6, eigen_quad_p310_2_7,
  eigen_quad_p310_2_8, eigen_quad_p310_2_9, eigen_quad_p310_2_10, eigen_quad_p310_3_2, eigen_quad_p310_3_3, eigen_quad_p310_3_4,
  eigen_quad_p310_3_5, eigen_quad_p310_3_6, eigen_quad_p310_3_7, eigen_quad_p310_3_8, eigen_quad_p310_3_9, eigen_quad_p310_3_10,

  eigen_quad_p42_2_2, eigen_quad_p42_3_2, eigen_quad_p42_4_2, eigen_quad_p43_2_2, eigen_quad_p43_2_3, eigen_quad_p43_3_2,
  eigen_quad_p43_3_3, eigen_quad_p43_4_2, eigen_quad_p43_4_3, eigen_quad_p44_2_2, eigen_quad_p44_2_3, eigen_quad_p44_2_4,
  eigen_quad_p44_3_2, eigen_quad_p44_3_3, eigen_quad_p44_3_4, eigen_quad_p44_4_2, eigen_quad_p44_4_3, eigen_quad_p44_4_4,
  eigen_quad_p45_2_2, eigen_quad_p45_2_3, eigen_quad_p45_2_4, eigen_quad_p45_2_5, eigen_quad_p45_3_2, eigen_quad_p45_3_3,
  eigen_quad_p45_3_4, eigen_quad_p45_3_5, eigen_quad_p45_4_2, eigen_quad_p45_4_3, eigen_quad_p45_4_4, eigen_quad_p45_4_5,
  eigen_quad_p46_2_2, eigen_quad_p46_2_3, eigen_quad_p46_2_4, eigen_quad_p46_2_5, eigen_quad_p46_2_6, eigen_quad_p46_3_2,
  eigen_quad_p46_3_3, eigen_quad_p46_3_4, eigen_quad_p46_3_5, eigen_quad_p46_3_6, eigen_quad_p46_4_2, eigen_quad_p46_4_3,
  eigen_quad_p46_4_4, eigen_quad_p46_4_5, eigen_quad_p46_4_6, eigen_quad_p47_2_2, eigen_quad_p47_2_3, eigen_quad_p47_2_4,
  eigen_quad_p47_2_5, eigen_quad_p47_2_6, eigen_quad_p47_2_7, eigen_quad_p47_3_2, eigen_quad_p47_3_3, eigen_quad_p47_3_4,
  eigen_quad_p47_3_5, eigen_quad_p47_3_6, eigen_quad_p47_3_7, eigen_quad_p47_4_2, eigen_quad_p47_4_3, eigen_quad_p47_4_4,
  eigen_quad_p47_4_5, eigen_quad_p47_4_6, eigen_quad_p47_4_7, eigen_quad_p48_2_2, eigen_quad_p48_2_3, eigen_quad_p48_2_4,
  eigen_quad_p48_2_5, eigen_quad_p48_2_6, eigen_quad_p48_2_7, eigen_quad_p48_2_8, eigen_quad_p48_3_2, eigen_quad_p48_3_3,
  eigen_quad_p48_3_4, eigen_quad_p48_3_5, eigen_quad_p48_3_6, eigen_quad_p48_3_7, eigen_quad_p48_3_8, eigen_quad_p48_4_2,
  eigen_quad_p48_4_3, eigen_quad_p48_4_4, eigen_quad_p48_4_5, eigen_quad_p48_4_6, eigen_quad_p48_4_7, eigen_quad_p48_4_8,
  eigen_quad_p49_2_2, eigen_quad_p49_2_3, eigen_quad_p49_2_4, eigen_quad_p49_2_5, eigen_quad_p49_2_6, eigen_quad_p49_2_7,
  eigen_quad_p49_2_8, eigen_quad_p49_2_9, eigen_quad_p49_3_2, eigen_quad_p49_3_3, eigen_quad_p49_3_4, eigen_quad_p49_3_5,
  eigen_quad_p49_3_6, eigen_quad_p49_3_7, eigen_quad_p49_3_8, eigen_quad_p49_3_9, eigen_quad_p49_4_2, eigen_quad_p49_4_3,
  eigen_quad_p49_4_4, eigen_quad_p49_4_5, eigen_quad_p49_4_6, eigen_quad_p49_4_7, eigen_quad_p49_4_8, eigen_quad_p49_4_9,
  eigen_quad_p410_2_2, eigen_quad_p410_2_3, eigen_quad_p410_2_4, eigen_quad_p410_2_5, eigen_quad_p410_2_6, eigen_quad_p410_2_7,
  eigen_quad_p410_2_8, eigen_quad_p410_2_9, eigen_quad_p410_2_10, eigen_quad_p410_3_2, eigen_quad_p410_3_3, eigen_quad_p410_3_4,
  eigen_quad_p410_3_5, eigen_quad_p410_3_6, eigen_quad_p410_3_7, eigen_quad_p410_3_8, eigen_quad_p410_3_9, eigen_quad_p410_3_10,
  eigen_quad_p410_4_2, eigen_quad_p410_4_3, eigen_quad_p410_4_4, eigen_quad_p410_4_5, eigen_quad_p410_4_6, eigen_quad_p410_4_7,
  eigen_quad_p410_4_8, eigen_quad_p410_4_9, eigen_quad_p410_4_10,

  eigen_quad_p52_2_2, eigen_quad_p52_3_2, eigen_quad_p52_4_2, eigen_quad_p52_5_2, eigen_quad_p53_2_2, eigen_quad_p53_2_3, eigen_quad_p53_3_2, eigen_quad_p53_3_3, eigen_quad_p53_4_2, eigen_quad_p53_4_3, eigen_quad_p53_5_2, eigen_quad_p53_5_3, eigen_quad_p54_2_2, eigen_quad_p54_2_3, eigen_quad_p54_2_4, eigen_quad_p54_3_2, eigen_quad_p54_3_3, eigen_quad_p54_3_4, eigen_quad_p54_4_2, eigen_quad_p54_4_3, eigen_quad_p54_4_4, eigen_quad_p54_5_2, eigen_quad_p54_5_3, eigen_quad_p54_5_4, eigen_quad_p55_2_2, eigen_quad_p55_2_3, eigen_quad_p55_2_4, eigen_quad_p55_2_5, eigen_quad_p55_3_2, eigen_quad_p55_3_3, eigen_quad_p55_3_4, eigen_quad_p55_3_5, eigen_quad_p55_4_2, eigen_quad_p55_4_3, eigen_quad_p55_4_4, eigen_quad_p55_4_5, eigen_quad_p55_5_2, eigen_quad_p55_5_3, eigen_quad_p55_5_4, eigen_quad_p55_5_5, eigen_quad_p56_2_2, eigen_quad_p56_2_3, eigen_quad_p56_2_4, eigen_quad_p56_2_5, eigen_quad_p56_2_6, eigen_quad_p56_3_2, eigen_quad_p56_3_3, eigen_quad_p56_3_4, eigen_quad_p56_3_5, eigen_quad_p56_3_6, eigen_quad_p56_4_2, eigen_quad_p56_4_3, eigen_quad_p56_4_4, eigen_quad_p56_4_5, eigen_quad_p56_4_6, eigen_quad_p56_5_2, eigen_quad_p56_5_3, eigen_quad_p56_5_4, eigen_quad_p56_5_5, eigen_quad_p56_5_6, eigen_quad_p57_2_2, eigen_quad_p57_2_3, eigen_quad_p57_2_4, eigen_quad_p57_2_5, eigen_quad_p57_2_6, eigen_quad_p57_2_7, eigen_quad_p57_3_2, eigen_quad_p57_3_3, eigen_quad_p57_3_4, eigen_quad_p57_3_5, eigen_quad_p57_3_6, eigen_quad_p57_3_7, eigen_quad_p57_4_2, eigen_quad_p57_4_3, eigen_quad_p57_4_4, eigen_quad_p57_4_5, eigen_quad_p57_4_6, eigen_quad_p57_4_7, eigen_quad_p57_5_2, eigen_quad_p57_5_3, eigen_quad_p57_5_4, eigen_quad_p57_5_5, eigen_quad_p57_5_6, eigen_quad_p57_5_7, eigen_quad_p58_2_2, eigen_quad_p58_2_3, eigen_quad_p58_2_4, eigen_quad_p58_2_5, eigen_quad_p58_2_6, eigen_quad_p58_2_7, eigen_quad_p58_2_8, eigen_quad_p58_3_2, eigen_quad_p58_3_3, eigen_quad_p58_3_4, eigen_quad_p58_3_5, eigen_quad_p58_3_6, eigen_quad_p58_3_7, eigen_quad_p58_3_8, eigen_quad_p58_4_2, eigen_quad_p58_4_3, eigen_quad_p58_4_4, eigen_quad_p58_4_5, eigen_quad_p58_4_6, eigen_quad_p58_4_7, eigen_quad_p58_4_8, eigen_quad_p58_5_2, eigen_quad_p58_5_3, eigen_quad_p58_5_4, eigen_quad_p58_5_5, eigen_quad_p58_5_6, eigen_quad_p58_5_7, eigen_quad_p58_5_8, eigen_quad_p59_2_2, eigen_quad_p59_2_3, eigen_quad_p59_2_4, eigen_quad_p59_2_5, eigen_quad_p59_2_6, eigen_quad_p59_2_7, eigen_quad_p59_2_8, eigen_quad_p59_2_9, eigen_quad_p59_3_2, eigen_quad_p59_3_3, eigen_quad_p59_3_4, eigen_quad_p59_3_5, eigen_quad_p59_3_6, eigen_quad_p59_3_7, eigen_quad_p59_3_8, eigen_quad_p59_3_9, eigen_quad_p59_4_2, eigen_quad_p59_4_3, eigen_quad_p59_4_4, eigen_quad_p59_4_5, eigen_quad_p59_4_6, eigen_quad_p59_4_7, eigen_quad_p59_4_8, eigen_quad_p59_4_9, eigen_quad_p59_5_2, eigen_quad_p59_5_3, eigen_quad_p59_5_4, eigen_quad_p59_5_5, eigen_quad_p59_5_6, eigen_quad_p59_5_7, eigen_quad_p59_5_8, eigen_quad_p59_5_9, eigen_quad_p510_2_2, eigen_quad_p510_2_3, eigen_quad_p510_2_4, eigen_quad_p510_2_5, eigen_quad_p510_2_6, eigen_quad_p510_2_7, eigen_quad_p510_2_8, eigen_quad_p510_2_9, eigen_quad_p510_2_10, eigen_quad_p510_3_2, eigen_quad_p510_3_3, eigen_quad_p510_3_4, eigen_quad_p510_3_5, eigen_quad_p510_3_6, eigen_quad_p510_3_7, eigen_quad_p510_3_8, eigen_quad_p510_3_9, eigen_quad_p510_3_10, eigen_quad_p510_4_2, eigen_quad_p510_4_3, eigen_quad_p510_4_4, eigen_quad_p510_4_5, eigen_quad_p510_4_6, eigen_quad_p510_4_7, eigen_quad_p510_4_8, eigen_quad_p510_4_9, eigen_quad_p510_4_10, eigen_quad_p510_5_2, eigen_quad_p510_5_3, eigen_quad_p510_5_4, eigen_quad_p510_5_5, eigen_quad_p510_5_6, eigen_quad_p510_5_7, eigen_quad_p510_5_8, eigen_quad_p510_5_9, eigen_quad_p510_5_10,

  eigen_quad_p62_2_2, eigen_quad_p62_3_2, eigen_quad_p62_4_2, eigen_quad_p62_5_2, eigen_quad_p62_6_2, eigen_quad_p63_2_2, eigen_quad_p63_2_3, eigen_quad_p63_3_2, eigen_quad_p63_3_3, eigen_quad_p63_4_2, eigen_quad_p63_4_3, eigen_quad_p63_5_2, eigen_quad_p63_5_3, eigen_quad_p63_6_2, eigen_quad_p63_6_3, eigen_quad_p64_2_2, eigen_quad_p64_2_3, eigen_quad_p64_2_4, eigen_quad_p64_3_2, eigen_quad_p64_3_3, eigen_quad_p64_3_4, eigen_quad_p64_4_2, eigen_quad_p64_4_3, eigen_quad_p64_4_4, eigen_quad_p64_5_2, eigen_quad_p64_5_3, eigen_quad_p64_5_4, eigen_quad_p64_6_2, eigen_quad_p64_6_3, eigen_quad_p64_6_4, eigen_quad_p65_2_2, eigen_quad_p65_2_3, eigen_quad_p65_2_4, eigen_quad_p65_2_5, eigen_quad_p65_3_2, eigen_quad_p65_3_3, eigen_quad_p65_3_4, eigen_quad_p65_3_5, eigen_quad_p65_4_2, eigen_quad_p65_4_3, eigen_quad_p65_4_4, eigen_quad_p65_4_5, eigen_quad_p65_5_2, eigen_quad_p65_5_3, eigen_quad_p65_5_4, eigen_quad_p65_5_5, eigen_quad_p65_6_2, eigen_quad_p65_6_3, eigen_quad_p65_6_4, eigen_quad_p65_6_5, eigen_quad_p66_2_2, eigen_quad_p66_2_3, eigen_quad_p66_2_4, eigen_quad_p66_2_5, eigen_quad_p66_2_6, eigen_quad_p66_3_2, eigen_quad_p66_3_3, eigen_quad_p66_3_4, eigen_quad_p66_3_5, eigen_quad_p66_3_6, eigen_quad_p66_4_2, eigen_quad_p66_4_3, eigen_quad_p66_4_4, eigen_quad_p66_4_5, eigen_quad_p66_4_6, eigen_quad_p66_5_2, eigen_quad_p66_5_3, eigen_quad_p66_5_4, eigen_quad_p66_5_5, eigen_quad_p66_5_6, eigen_quad_p66_6_2, eigen_quad_p66_6_3, eigen_quad_p66_6_4, eigen_quad_p66_6_5, eigen_quad_p66_6_6, eigen_quad_p67_2_2, eigen_quad_p67_2_3, eigen_quad_p67_2_4, eigen_quad_p67_2_5, eigen_quad_p67_2_6, eigen_quad_p67_2_7, eigen_quad_p67_3_2, eigen_quad_p67_3_3, eigen_quad_p67_3_4, eigen_quad_p67_3_5, eigen_quad_p67_3_6, eigen_quad_p67_3_7, eigen_quad_p67_4_2, eigen_quad_p67_4_3, eigen_quad_p67_4_4, eigen_quad_p67_4_5, eigen_quad_p67_4_6, eigen_quad_p67_4_7, eigen_quad_p67_5_2, eigen_quad_p67_5_3, eigen_quad_p67_5_4, eigen_quad_p67_5_5, eigen_quad_p67_5_6, eigen_quad_p67_5_7, eigen_quad_p67_6_2, eigen_quad_p67_6_3, eigen_quad_p67_6_4, eigen_quad_p67_6_5, eigen_quad_p67_6_6, eigen_quad_p67_6_7, eigen_quad_p68_2_2, eigen_quad_p68_2_3, eigen_quad_p68_2_4, eigen_quad_p68_2_5, eigen_quad_p68_2_6, eigen_quad_p68_2_7, eigen_quad_p68_2_8, eigen_quad_p68_3_2, eigen_quad_p68_3_3, eigen_quad_p68_3_4, eigen_quad_p68_3_5, eigen_quad_p68_3_6, eigen_quad_p68_3_7, eigen_quad_p68_3_8, eigen_quad_p68_4_2, eigen_quad_p68_4_3, eigen_quad_p68_4_4, eigen_quad_p68_4_5, eigen_quad_p68_4_6, eigen_quad_p68_4_7, eigen_quad_p68_4_8, eigen_quad_p68_5_2, eigen_quad_p68_5_3, eigen_quad_p68_5_4, eigen_quad_p68_5_5, eigen_quad_p68_5_6, eigen_quad_p68_5_7, eigen_quad_p68_5_8, eigen_quad_p68_6_2, eigen_quad_p68_6_3, eigen_quad_p68_6_4, eigen_quad_p68_6_5, eigen_quad_p68_6_6, eigen_quad_p68_6_7, eigen_quad_p68_6_8, eigen_quad_p69_2_2, eigen_quad_p69_2_3, eigen_quad_p69_2_4, eigen_quad_p69_2_5, eigen_quad_p69_2_6, eigen_quad_p69_2_7, eigen_quad_p69_2_8, eigen_quad_p69_2_9, eigen_quad_p69_3_2, eigen_quad_p69_3_3, eigen_quad_p69_3_4, eigen_quad_p69_3_5, eigen_quad_p69_3_6, eigen_quad_p69_3_7, eigen_quad_p69_3_8, eigen_quad_p69_3_9, eigen_quad_p69_4_2, eigen_quad_p69_4_3, eigen_quad_p69_4_4, eigen_quad_p69_4_5, eigen_quad_p69_4_6, eigen_quad_p69_4_7, eigen_quad_p69_4_8, eigen_quad_p69_4_9, eigen_quad_p69_5_2, eigen_quad_p69_5_3, eigen_quad_p69_5_4, eigen_quad_p69_5_5, eigen_quad_p69_5_6, eigen_quad_p69_5_7, eigen_quad_p69_5_8, eigen_quad_p69_5_9, eigen_quad_p69_6_2, eigen_quad_p69_6_3, eigen_quad_p69_6_4, eigen_quad_p69_6_5, eigen_quad_p69_6_6, eigen_quad_p69_6_7, eigen_quad_p69_6_8, eigen_quad_p69_6_9, eigen_quad_p610_2_2, eigen_quad_p610_2_3, eigen_quad_p610_2_4, eigen_quad_p610_2_5, eigen_quad_p610_2_6, eigen_quad_p610_2_7, eigen_quad_p610_2_8, eigen_quad_p610_2_9, eigen_quad_p610_2_10, eigen_quad_p610_3_2, eigen_quad_p610_3_3, eigen_quad_p610_3_4, eigen_quad_p610_3_5, eigen_quad_p610_3_6, eigen_quad_p610_3_7, eigen_quad_p610_3_8, eigen_quad_p610_3_9, eigen_quad_p610_3_10, eigen_quad_p610_4_2, eigen_quad_p610_4_3, eigen_quad_p610_4_4, eigen_quad_p610_4_5, eigen_quad_p610_4_6, eigen_quad_p610_4_7, eigen_quad_p610_4_8, eigen_quad_p610_4_9, eigen_quad_p610_4_10, eigen_quad_p610_5_2, eigen_quad_p610_5_3, eigen_quad_p610_5_4, eigen_quad_p610_5_5, eigen_quad_p610_5_6, eigen_quad_p610_5_7, eigen_quad_p610_5_8, eigen_quad_p610_5_9, eigen_quad_p610_5_10, eigen_quad_p610_6_2, eigen_quad_p610_6_3, eigen_quad_p610_6_4, eigen_quad_p610_6_5, eigen_quad_p610_6_6, eigen_quad_p610_6_7, eigen_quad_p610_6_8, eigen_quad_p610_6_9, eigen_quad_p610_6_10,

  eigen_quad_p72_2_2, eigen_quad_p72_3_2, eigen_quad_p72_4_2, eigen_quad_p72_5_2, eigen_quad_p72_6_2, eigen_quad_p72_7_2, eigen_quad_p73_2_2, eigen_quad_p73_2_3, eigen_quad_p73_3_2, eigen_quad_p73_3_3, eigen_quad_p73_4_2, eigen_quad_p73_4_3, eigen_quad_p73_5_2, eigen_quad_p73_5_3, eigen_quad_p73_6_2, eigen_quad_p73_6_3, eigen_quad_p73_7_2, eigen_quad_p73_7_3, eigen_quad_p74_2_2, eigen_quad_p74_2_3, eigen_quad_p74_2_4, eigen_quad_p74_3_2, eigen_quad_p74_3_3, eigen_quad_p74_3_4, eigen_quad_p74_4_2, eigen_quad_p74_4_3, eigen_quad_p74_4_4, eigen_quad_p74_5_2, eigen_quad_p74_5_3, eigen_quad_p74_5_4, eigen_quad_p74_6_2, eigen_quad_p74_6_3, eigen_quad_p74_6_4, eigen_quad_p74_7_2, eigen_quad_p74_7_3, eigen_quad_p74_7_4, eigen_quad_p75_2_2, eigen_quad_p75_2_3, eigen_quad_p75_2_4, eigen_quad_p75_2_5, eigen_quad_p75_3_2, eigen_quad_p75_3_3, eigen_quad_p75_3_4, eigen_quad_p75_3_5, eigen_quad_p75_4_2, eigen_quad_p75_4_3, eigen_quad_p75_4_4, eigen_quad_p75_4_5, eigen_quad_p75_5_2, eigen_quad_p75_5_3, eigen_quad_p75_5_4, eigen_quad_p75_5_5, eigen_quad_p75_6_2, eigen_quad_p75_6_3, eigen_quad_p75_6_4, eigen_quad_p75_6_5, eigen_quad_p75_7_2, eigen_quad_p75_7_3, eigen_quad_p75_7_4, eigen_quad_p75_7_5, eigen_quad_p76_2_2, eigen_quad_p76_2_3, eigen_quad_p76_2_4, eigen_quad_p76_2_5, eigen_quad_p76_2_6, eigen_quad_p76_3_2, eigen_quad_p76_3_3, eigen_quad_p76_3_4, eigen_quad_p76_3_5, eigen_quad_p76_3_6, eigen_quad_p76_4_2, eigen_quad_p76_4_3, eigen_quad_p76_4_4, eigen_quad_p76_4_5, eigen_quad_p76_4_6, eigen_quad_p76_5_2, eigen_quad_p76_5_3, eigen_quad_p76_5_4, eigen_quad_p76_5_5, eigen_quad_p76_5_6, eigen_quad_p76_6_2, eigen_quad_p76_6_3, eigen_quad_p76_6_4, eigen_quad_p76_6_5, eigen_quad_p76_6_6, eigen_quad_p76_7_2, eigen_quad_p76_7_3, eigen_quad_p76_7_4, eigen_quad_p76_7_5, eigen_quad_p76_7_6, eigen_quad_p77_2_2, eigen_quad_p77_2_3, eigen_quad_p77_2_4, eigen_quad_p77_2_5, eigen_quad_p77_2_6, eigen_quad_p77_2_7, eigen_quad_p77_3_2, eigen_quad_p77_3_3, eigen_quad_p77_3_4, eigen_quad_p77_3_5, eigen_quad_p77_3_6, eigen_quad_p77_3_7, eigen_quad_p77_4_2, eigen_quad_p77_4_3, eigen_quad_p77_4_4, eigen_quad_p77_4_5, eigen_quad_p77_4_6, eigen_quad_p77_4_7, eigen_quad_p77_5_2, eigen_quad_p77_5_3, eigen_quad_p77_5_4, eigen_quad_p77_5_5, eigen_quad_p77_5_6, eigen_quad_p77_5_7, eigen_quad_p77_6_2, eigen_quad_p77_6_3, eigen_quad_p77_6_4, eigen_quad_p77_6_5, eigen_quad_p77_6_6, eigen_quad_p77_6_7, eigen_quad_p77_7_2, eigen_quad_p77_7_3, eigen_quad_p77_7_4, eigen_quad_p77_7_5, eigen_quad_p77_7_6, eigen_quad_p77_7_7, eigen_quad_p78_2_2, eigen_quad_p78_2_3, eigen_quad_p78_2_4, eigen_quad_p78_2_5, eigen_quad_p78_2_6, eigen_quad_p78_2_7, eigen_quad_p78_2_8, eigen_quad_p78_3_2, eigen_quad_p78_3_3, eigen_quad_p78_3_4, eigen_quad_p78_3_5, eigen_quad_p78_3_6, eigen_quad_p78_3_7, eigen_quad_p78_3_8, eigen_quad_p78_4_2, eigen_quad_p78_4_3, eigen_quad_p78_4_4, eigen_quad_p78_4_5, eigen_quad_p78_4_6, eigen_quad_p78_4_7, eigen_quad_p78_4_8, eigen_quad_p78_5_2, eigen_quad_p78_5_3, eigen_quad_p78_5_4, eigen_quad_p78_5_5, eigen_quad_p78_5_6, eigen_quad_p78_5_7, eigen_quad_p78_5_8, eigen_quad_p78_6_2, eigen_quad_p78_6_3, eigen_quad_p78_6_4, eigen_quad_p78_6_5, eigen_quad_p78_6_6, eigen_quad_p78_6_7, eigen_quad_p78_6_8, eigen_quad_p78_7_2, eigen_quad_p78_7_3, eigen_quad_p78_7_4, eigen_quad_p78_7_5, eigen_quad_p78_7_6, eigen_quad_p78_7_7, eigen_quad_p78_7_8, eigen_quad_p79_2_2, eigen_quad_p79_2_3, eigen_quad_p79_2_4, eigen_quad_p79_2_5, eigen_quad_p79_2_6, eigen_quad_p79_2_7, eigen_quad_p79_2_8, eigen_quad_p79_2_9, eigen_quad_p79_3_2, eigen_quad_p79_3_3, eigen_quad_p79_3_4, eigen_quad_p79_3_5, eigen_quad_p79_3_6, eigen_quad_p79_3_7, eigen_quad_p79_3_8, eigen_quad_p79_3_9, eigen_quad_p79_4_2, eigen_quad_p79_4_3, eigen_quad_p79_4_4, eigen_quad_p79_4_5, eigen_quad_p79_4_6, eigen_quad_p79_4_7, eigen_quad_p79_4_8, eigen_quad_p79_4_9, eigen_quad_p79_5_2, eigen_quad_p79_5_3, eigen_quad_p79_5_4, eigen_quad_p79_5_5, eigen_quad_p79_5_6, eigen_quad_p79_5_7, eigen_quad_p79_5_8, eigen_quad_p79_5_9, eigen_quad_p79_6_2, eigen_quad_p79_6_3, eigen_quad_p79_6_4, eigen_quad_p79_6_5, eigen_quad_p79_6_6, eigen_quad_p79_6_7, eigen_quad_p79_6_8, eigen_quad_p79_6_9, eigen_quad_p79_7_2, eigen_quad_p79_7_3, eigen_quad_p79_7_4, eigen_quad_p79_7_5, eigen_quad_p79_7_6, eigen_quad_p79_7_7, eigen_quad_p79_7_8, eigen_quad_p79_7_9, eigen_quad_p710_2_2, eigen_quad_p710_2_3, eigen_quad_p710_2_4, eigen_quad_p710_2_5, eigen_quad_p710_2_6, eigen_quad_p710_2_7, eigen_quad_p710_2_8, eigen_quad_p710_2_9, eigen_quad_p710_2_10, eigen_quad_p710_3_2, eigen_quad_p710_3_3, eigen_quad_p710_3_4, eigen_quad_p710_3_5, eigen_quad_p710_3_6, eigen_quad_p710_3_7, eigen_quad_p710_3_8, eigen_quad_p710_3_9, eigen_quad_p710_3_10, eigen_quad_p710_4_2, eigen_quad_p710_4_3, eigen_quad_p710_4_4, eigen_quad_p710_4_5, eigen_quad_p710_4_6, eigen_quad_p710_4_7, eigen_quad_p710_4_8, eigen_quad_p710_4_9, eigen_quad_p710_4_10, eigen_quad_p710_5_2, eigen_quad_p710_5_3, eigen_quad_p710_5_4, eigen_quad_p710_5_5, eigen_quad_p710_5_6, eigen_quad_p710_5_7, eigen_quad_p710_5_8, eigen_quad_p710_5_9, eigen_quad_p710_5_10, eigen_quad_p710_6_2, eigen_quad_p710_6_3, eigen_quad_p710_6_4, eigen_quad_p710_6_5, eigen_quad_p710_6_6, eigen_quad_p710_6_7, eigen_quad_p710_6_8, eigen_quad_p710_6_9, eigen_quad_p710_6_10, eigen_quad_p710_7_2, eigen_quad_p710_7_3, eigen_quad_p710_7_4, eigen_quad_p710_7_5, eigen_quad_p710_7_6, eigen_quad_p710_7_7, eigen_quad_p710_7_8, eigen_quad_p710_7_9, eigen_quad_p710_7_10,

  eigen_quad_p82_2_2, eigen_quad_p82_3_2, eigen_quad_p82_4_2, eigen_quad_p82_5_2, eigen_quad_p82_6_2, eigen_quad_p82_7_2, eigen_quad_p82_8_2, eigen_quad_p83_2_2, eigen_quad_p83_2_3, eigen_quad_p83_3_2, eigen_quad_p83_3_3, eigen_quad_p83_4_2, eigen_quad_p83_4_3, eigen_quad_p83_5_2, eigen_quad_p83_5_3, eigen_quad_p83_6_2, eigen_quad_p83_6_3, eigen_quad_p83_7_2, eigen_quad_p83_7_3, eigen_quad_p83_8_2, eigen_quad_p83_8_3, eigen_quad_p84_2_2, eigen_quad_p84_2_3, eigen_quad_p84_2_4, eigen_quad_p84_3_2, eigen_quad_p84_3_3, eigen_quad_p84_3_4, eigen_quad_p84_4_2, eigen_quad_p84_4_3, eigen_quad_p84_4_4, eigen_quad_p84_5_2, eigen_quad_p84_5_3, eigen_quad_p84_5_4, eigen_quad_p84_6_2, eigen_quad_p84_6_3, eigen_quad_p84_6_4, eigen_quad_p84_7_2, eigen_quad_p84_7_3, eigen_quad_p84_7_4, eigen_quad_p84_8_2, eigen_quad_p84_8_3, eigen_quad_p84_8_4, eigen_quad_p85_2_2, eigen_quad_p85_2_3, eigen_quad_p85_2_4, eigen_quad_p85_2_5, eigen_quad_p85_3_2, eigen_quad_p85_3_3, eigen_quad_p85_3_4, eigen_quad_p85_3_5, eigen_quad_p85_4_2, eigen_quad_p85_4_3, eigen_quad_p85_4_4, eigen_quad_p85_4_5, eigen_quad_p85_5_2, eigen_quad_p85_5_3, eigen_quad_p85_5_4, eigen_quad_p85_5_5, eigen_quad_p85_6_2, eigen_quad_p85_6_3, eigen_quad_p85_6_4, eigen_quad_p85_6_5, eigen_quad_p85_7_2, eigen_quad_p85_7_3, eigen_quad_p85_7_4, eigen_quad_p85_7_5, eigen_quad_p85_8_2, eigen_quad_p85_8_3, eigen_quad_p85_8_4, eigen_quad_p85_8_5, eigen_quad_p86_2_2, eigen_quad_p86_2_3, eigen_quad_p86_2_4, eigen_quad_p86_2_5, eigen_quad_p86_2_6, eigen_quad_p86_3_2, eigen_quad_p86_3_3, eigen_quad_p86_3_4, eigen_quad_p86_3_5, eigen_quad_p86_3_6, eigen_quad_p86_4_2, eigen_quad_p86_4_3, eigen_quad_p86_4_4, eigen_quad_p86_4_5, eigen_quad_p86_4_6, eigen_quad_p86_5_2, eigen_quad_p86_5_3, eigen_quad_p86_5_4, eigen_quad_p86_5_5, eigen_quad_p86_5_6, eigen_quad_p86_6_2, eigen_quad_p86_6_3, eigen_quad_p86_6_4, eigen_quad_p86_6_5, eigen_quad_p86_6_6, eigen_quad_p86_7_2, eigen_quad_p86_7_3, eigen_quad_p86_7_4, eigen_quad_p86_7_5, eigen_quad_p86_7_6, eigen_quad_p86_8_2, eigen_quad_p86_8_3, eigen_quad_p86_8_4, eigen_quad_p86_8_5, eigen_quad_p86_8_6, eigen_quad_p87_2_2, eigen_quad_p87_2_3, eigen_quad_p87_2_4, eigen_quad_p87_2_5, eigen_quad_p87_2_6, eigen_quad_p87_2_7, eigen_quad_p87_3_2, eigen_quad_p87_3_3, eigen_quad_p87_3_4, eigen_quad_p87_3_5, eigen_quad_p87_3_6, eigen_quad_p87_3_7, eigen_quad_p87_4_2, eigen_quad_p87_4_3, eigen_quad_p87_4_4, eigen_quad_p87_4_5, eigen_quad_p87_4_6, eigen_quad_p87_4_7, eigen_quad_p87_5_2, eigen_quad_p87_5_3, eigen_quad_p87_5_4, eigen_quad_p87_5_5, eigen_quad_p87_5_6, eigen_quad_p87_5_7, eigen_quad_p87_6_2, eigen_quad_p87_6_3, eigen_quad_p87_6_4, eigen_quad_p87_6_5, eigen_quad_p87_6_6, eigen_quad_p87_6_7, eigen_quad_p87_7_2, eigen_quad_p87_7_3, eigen_quad_p87_7_4, eigen_quad_p87_7_5, eigen_quad_p87_7_6, eigen_quad_p87_7_7, eigen_quad_p87_8_2, eigen_quad_p87_8_3, eigen_quad_p87_8_4, eigen_quad_p87_8_5, eigen_quad_p87_8_6, eigen_quad_p87_8_7, eigen_quad_p88_2_2, eigen_quad_p88_2_3, eigen_quad_p88_2_4, eigen_quad_p88_2_5, eigen_quad_p88_2_6, eigen_quad_p88_2_7, eigen_quad_p88_2_8, eigen_quad_p88_3_2, eigen_quad_p88_3_3, eigen_quad_p88_3_4, eigen_quad_p88_3_5, eigen_quad_p88_3_6, eigen_quad_p88_3_7, eigen_quad_p88_3_8, eigen_quad_p88_4_2, eigen_quad_p88_4_3, eigen_quad_p88_4_4, eigen_quad_p88_4_5, eigen_quad_p88_4_6, eigen_quad_p88_4_7, eigen_quad_p88_4_8, eigen_quad_p88_5_2, eigen_quad_p88_5_3, eigen_quad_p88_5_4, eigen_quad_p88_5_5, eigen_quad_p88_5_6, eigen_quad_p88_5_7, eigen_quad_p88_5_8, eigen_quad_p88_6_2, eigen_quad_p88_6_3, eigen_quad_p88_6_4, eigen_quad_p88_6_5, eigen_quad_p88_6_6, eigen_quad_p88_6_7, eigen_quad_p88_6_8, eigen_quad_p88_7_2, eigen_quad_p88_7_3, eigen_quad_p88_7_4, eigen_quad_p88_7_5, eigen_quad_p88_7_6, eigen_quad_p88_7_7, eigen_quad_p88_7_8, eigen_quad_p88_8_2, eigen_quad_p88_8_3, eigen_quad_p88_8_4, eigen_quad_p88_8_5, eigen_quad_p88_8_6, eigen_quad_p88_8_7, eigen_quad_p88_8_8, eigen_quad_p89_2_2, eigen_quad_p89_2_3, eigen_quad_p89_2_4, eigen_quad_p89_2_5, eigen_quad_p89_2_6, eigen_quad_p89_2_7, eigen_quad_p89_2_8, eigen_quad_p89_2_9, eigen_quad_p89_3_2, eigen_quad_p89_3_3, eigen_quad_p89_3_4, eigen_quad_p89_3_5, eigen_quad_p89_3_6, eigen_quad_p89_3_7, eigen_quad_p89_3_8, eigen_quad_p89_3_9, eigen_quad_p89_4_2, eigen_quad_p89_4_3, eigen_quad_p89_4_4, eigen_quad_p89_4_5, eigen_quad_p89_4_6, eigen_quad_p89_4_7, eigen_quad_p89_4_8, eigen_quad_p89_4_9, eigen_quad_p89_5_2, eigen_quad_p89_5_3, eigen_quad_p89_5_4, eigen_quad_p89_5_5, eigen_quad_p89_5_6, eigen_quad_p89_5_7, eigen_quad_p89_5_8, eigen_quad_p89_5_9, eigen_quad_p89_6_2, eigen_quad_p89_6_3, eigen_quad_p89_6_4, eigen_quad_p89_6_5, eigen_quad_p89_6_6, eigen_quad_p89_6_7, eigen_quad_p89_6_8, eigen_quad_p89_6_9, eigen_quad_p89_7_2, eigen_quad_p89_7_3, eigen_quad_p89_7_4, eigen_quad_p89_7_5, eigen_quad_p89_7_6, eigen_quad_p89_7_7, eigen_quad_p89_7_8, eigen_quad_p89_7_9, eigen_quad_p89_8_2, eigen_quad_p89_8_3, eigen_quad_p89_8_4, eigen_quad_p89_8_5, eigen_quad_p89_8_6, eigen_quad_p89_8_7, eigen_quad_p89_8_8, eigen_quad_p89_8_9, eigen_quad_p810_2_2, eigen_quad_p810_2_3, eigen_quad_p810_2_4, eigen_quad_p810_2_5, eigen_quad_p810_2_6, eigen_quad_p810_2_7, eigen_quad_p810_2_8, eigen_quad_p810_2_9, eigen_quad_p810_2_10, eigen_quad_p810_3_2, eigen_quad_p810_3_3, eigen_quad_p810_3_4, eigen_quad_p810_3_5, eigen_quad_p810_3_6, eigen_quad_p810_3_7, eigen_quad_p810_3_8, eigen_quad_p810_3_9, eigen_quad_p810_3_10, eigen_quad_p810_4_2, eigen_quad_p810_4_3, eigen_quad_p810_4_4, eigen_quad_p810_4_5, eigen_quad_p810_4_6, eigen_quad_p810_4_7, eigen_quad_p810_4_8, eigen_quad_p810_4_9, eigen_quad_p810_4_10, eigen_quad_p810_5_2, eigen_quad_p810_5_3, eigen_quad_p810_5_4, eigen_quad_p810_5_5, eigen_quad_p810_5_6, eigen_quad_p810_5_7, eigen_quad_p810_5_8, eigen_quad_p810_5_9, eigen_quad_p810_5_10, eigen_quad_p810_6_2, eigen_quad_p810_6_3, eigen_quad_p810_6_4, eigen_quad_p810_6_5, eigen_quad_p810_6_6, eigen_quad_p810_6_7, eigen_quad_p810_6_8, eigen_quad_p810_6_9, eigen_quad_p810_6_10, eigen_quad_p810_7_2, eigen_quad_p810_7_3, eigen_quad_p810_7_4, eigen_quad_p810_7_5, eigen_quad_p810_7_6, eigen_quad_p810_7_7, eigen_quad_p810_7_8, eigen_quad_p810_7_9, eigen_quad_p810_7_10, eigen_quad_p810_8_2, eigen_quad_p810_8_3, eigen_quad_p810_8_4, eigen_quad_p810_8_5, eigen_quad_p810_8_6, eigen_quad_p810_8_7, eigen_quad_p810_8_8, eigen_quad_p810_8_9, eigen_quad_p810_8_10,

  eigen_quad_p92_2_2, eigen_quad_p92_3_2, eigen_quad_p92_4_2, eigen_quad_p92_5_2, eigen_quad_p92_6_2, eigen_quad_p92_7_2, eigen_quad_p92_8_2, eigen_quad_p92_9_2, eigen_quad_p93_2_2, eigen_quad_p93_2_3, eigen_quad_p93_3_2, eigen_quad_p93_3_3, eigen_quad_p93_4_2, eigen_quad_p93_4_3, eigen_quad_p93_5_2, eigen_quad_p93_5_3, eigen_quad_p93_6_2, eigen_quad_p93_6_3, eigen_quad_p93_7_2, eigen_quad_p93_7_3, eigen_quad_p93_8_2, eigen_quad_p93_8_3, eigen_quad_p93_9_2, eigen_quad_p93_9_3, eigen_quad_p94_2_2, eigen_quad_p94_2_3, eigen_quad_p94_2_4, eigen_quad_p94_3_2, eigen_quad_p94_3_3, eigen_quad_p94_3_4, eigen_quad_p94_4_2, eigen_quad_p94_4_3, eigen_quad_p94_4_4, eigen_quad_p94_5_2, eigen_quad_p94_5_3, eigen_quad_p94_5_4, eigen_quad_p94_6_2, eigen_quad_p94_6_3, eigen_quad_p94_6_4, eigen_quad_p94_7_2, eigen_quad_p94_7_3, eigen_quad_p94_7_4, eigen_quad_p94_8_2, eigen_quad_p94_8_3, eigen_quad_p94_8_4, eigen_quad_p94_9_2, eigen_quad_p94_9_3, eigen_quad_p94_9_4, eigen_quad_p95_2_2, eigen_quad_p95_2_3, eigen_quad_p95_2_4, eigen_quad_p95_2_5, eigen_quad_p95_3_2, eigen_quad_p95_3_3, eigen_quad_p95_3_4, eigen_quad_p95_3_5, eigen_quad_p95_4_2, eigen_quad_p95_4_3, eigen_quad_p95_4_4, eigen_quad_p95_4_5, eigen_quad_p95_5_2, eigen_quad_p95_5_3, eigen_quad_p95_5_4, eigen_quad_p95_5_5, eigen_quad_p95_6_2, eigen_quad_p95_6_3, eigen_quad_p95_6_4, eigen_quad_p95_6_5, eigen_quad_p95_7_2, eigen_quad_p95_7_3, eigen_quad_p95_7_4, eigen_quad_p95_7_5, eigen_quad_p95_8_2, eigen_quad_p95_8_3, eigen_quad_p95_8_4, eigen_quad_p95_8_5, eigen_quad_p95_9_2, eigen_quad_p95_9_3, eigen_quad_p95_9_4, eigen_quad_p95_9_5, eigen_quad_p96_2_2, eigen_quad_p96_2_3, eigen_quad_p96_2_4, eigen_quad_p96_2_5, eigen_quad_p96_2_6, eigen_quad_p96_3_2, eigen_quad_p96_3_3, eigen_quad_p96_3_4, eigen_quad_p96_3_5, eigen_quad_p96_3_6, eigen_quad_p96_4_2, eigen_quad_p96_4_3, eigen_quad_p96_4_4, eigen_quad_p96_4_5, eigen_quad_p96_4_6, eigen_quad_p96_5_2, eigen_quad_p96_5_3, eigen_quad_p96_5_4, eigen_quad_p96_5_5, eigen_quad_p96_5_6, eigen_quad_p96_6_2, eigen_quad_p96_6_3, eigen_quad_p96_6_4, eigen_quad_p96_6_5, eigen_quad_p96_6_6, eigen_quad_p96_7_2, eigen_quad_p96_7_3, eigen_quad_p96_7_4, eigen_quad_p96_7_5, eigen_quad_p96_7_6, eigen_quad_p96_8_2, eigen_quad_p96_8_3, eigen_quad_p96_8_4, eigen_quad_p96_8_5, eigen_quad_p96_8_6, eigen_quad_p96_9_2, eigen_quad_p96_9_3, eigen_quad_p96_9_4, eigen_quad_p96_9_5, eigen_quad_p96_9_6, eigen_quad_p97_2_2, eigen_quad_p97_2_3, eigen_quad_p97_2_4, eigen_quad_p97_2_5, eigen_quad_p97_2_6, eigen_quad_p97_2_7, eigen_quad_p97_3_2, eigen_quad_p97_3_3, eigen_quad_p97_3_4, eigen_quad_p97_3_5, eigen_quad_p97_3_6, eigen_quad_p97_3_7, eigen_quad_p97_4_2, eigen_quad_p97_4_3, eigen_quad_p97_4_4, eigen_quad_p97_4_5, eigen_quad_p97_4_6, eigen_quad_p97_4_7, eigen_quad_p97_5_2, eigen_quad_p97_5_3, eigen_quad_p97_5_4, eigen_quad_p97_5_5, eigen_quad_p97_5_6, eigen_quad_p97_5_7, eigen_quad_p97_6_2, eigen_quad_p97_6_3, eigen_quad_p97_6_4, eigen_quad_p97_6_5, eigen_quad_p97_6_6, eigen_quad_p97_6_7, eigen_quad_p97_7_2, eigen_quad_p97_7_3, eigen_quad_p97_7_4, eigen_quad_p97_7_5, eigen_quad_p97_7_6, eigen_quad_p97_7_7, eigen_quad_p97_8_2, eigen_quad_p97_8_3, eigen_quad_p97_8_4, eigen_quad_p97_8_5, eigen_quad_p97_8_6, eigen_quad_p97_8_7, eigen_quad_p97_9_2, eigen_quad_p97_9_3, eigen_quad_p97_9_4, eigen_quad_p97_9_5, eigen_quad_p97_9_6, eigen_quad_p97_9_7, eigen_quad_p98_2_2, eigen_quad_p98_2_3, eigen_quad_p98_2_4, eigen_quad_p98_2_5, eigen_quad_p98_2_6, eigen_quad_p98_2_7, eigen_quad_p98_2_8, eigen_quad_p98_3_2, eigen_quad_p98_3_3, eigen_quad_p98_3_4, eigen_quad_p98_3_5, eigen_quad_p98_3_6, eigen_quad_p98_3_7, eigen_quad_p98_3_8, eigen_quad_p98_4_2, eigen_quad_p98_4_3, eigen_quad_p98_4_4, eigen_quad_p98_4_5, eigen_quad_p98_4_6, eigen_quad_p98_4_7, eigen_quad_p98_4_8, eigen_quad_p98_5_2, eigen_quad_p98_5_3, eigen_quad_p98_5_4, eigen_quad_p98_5_5, eigen_quad_p98_5_6, eigen_quad_p98_5_7, eigen_quad_p98_5_8, eigen_quad_p98_6_2, eigen_quad_p98_6_3, eigen_quad_p98_6_4, eigen_quad_p98_6_5, eigen_quad_p98_6_6, eigen_quad_p98_6_7, eigen_quad_p98_6_8, eigen_quad_p98_7_2, eigen_quad_p98_7_3, eigen_quad_p98_7_4, eigen_quad_p98_7_5, eigen_quad_p98_7_6, eigen_quad_p98_7_7, eigen_quad_p98_7_8, eigen_quad_p98_8_2, eigen_quad_p98_8_3, eigen_quad_p98_8_4, eigen_quad_p98_8_5, eigen_quad_p98_8_6, eigen_quad_p98_8_7, eigen_quad_p98_8_8, eigen_quad_p98_9_2, eigen_quad_p98_9_3, eigen_quad_p98_9_4, eigen_quad_p98_9_5, eigen_quad_p98_9_6, eigen_quad_p98_9_7, eigen_quad_p98_9_8, eigen_quad_p99_2_2, eigen_quad_p99_2_3, eigen_quad_p99_2_4, eigen_quad_p99_2_5, eigen_quad_p99_2_6, eigen_quad_p99_2_7, eigen_quad_p99_2_8, eigen_quad_p99_2_9, eigen_quad_p99_3_2, eigen_quad_p99_3_3, eigen_quad_p99_3_4, eigen_quad_p99_3_5, eigen_quad_p99_3_6, eigen_quad_p99_3_7, eigen_quad_p99_3_8, eigen_quad_p99_3_9, eigen_quad_p99_4_2, eigen_quad_p99_4_3, eigen_quad_p99_4_4, eigen_quad_p99_4_5, eigen_quad_p99_4_6, eigen_quad_p99_4_7, eigen_quad_p99_4_8, eigen_quad_p99_4_9, eigen_quad_p99_5_2, eigen_quad_p99_5_3, eigen_quad_p99_5_4, eigen_quad_p99_5_5, eigen_quad_p99_5_6, eigen_quad_p99_5_7, eigen_quad_p99_5_8, eigen_quad_p99_5_9, eigen_quad_p99_6_2, eigen_quad_p99_6_3, eigen_quad_p99_6_4, eigen_quad_p99_6_5, eigen_quad_p99_6_6, eigen_quad_p99_6_7, eigen_quad_p99_6_8, eigen_quad_p99_6_9, eigen_quad_p99_7_2, eigen_quad_p99_7_3, eigen_quad_p99_7_4, eigen_quad_p99_7_5, eigen_quad_p99_7_6, eigen_quad_p99_7_7, eigen_quad_p99_7_8, eigen_quad_p99_7_9, eigen_quad_p99_8_2, eigen_quad_p99_8_3, eigen_quad_p99_8_4, eigen_quad_p99_8_5, eigen_quad_p99_8_6, eigen_quad_p99_8_7, eigen_quad_p99_8_8, eigen_quad_p99_8_9, eigen_quad_p99_9_2, eigen_quad_p99_9_3, eigen_quad_p99_9_4, eigen_quad_p99_9_5, eigen_quad_p99_9_6, eigen_quad_p99_9_7, eigen_quad_p99_9_8, eigen_quad_p99_9_9, eigen_quad_p910_2_2, eigen_quad_p910_2_3, eigen_quad_p910_2_4, eigen_quad_p910_2_5, eigen_quad_p910_2_6, eigen_quad_p910_2_7, eigen_quad_p910_2_8, eigen_quad_p910_2_9, eigen_quad_p910_2_10, eigen_quad_p910_3_2, eigen_quad_p910_3_3, eigen_quad_p910_3_4, eigen_quad_p910_3_5, eigen_quad_p910_3_6, eigen_quad_p910_3_7, eigen_quad_p910_3_8, eigen_quad_p910_3_9, eigen_quad_p910_3_10, eigen_quad_p910_4_2, eigen_quad_p910_4_3, eigen_quad_p910_4_4, eigen_quad_p910_4_5, eigen_quad_p910_4_6, eigen_quad_p910_4_7, eigen_quad_p910_4_8, eigen_quad_p910_4_9, eigen_quad_p910_4_10, eigen_quad_p910_5_2, eigen_quad_p910_5_3, eigen_quad_p910_5_4, eigen_quad_p910_5_5, eigen_quad_p910_5_6, eigen_quad_p910_5_7, eigen_quad_p910_5_8, eigen_quad_p910_5_9, eigen_quad_p910_5_10, eigen_quad_p910_6_2, eigen_quad_p910_6_3, eigen_quad_p910_6_4, eigen_quad_p910_6_5, eigen_quad_p910_6_6, eigen_quad_p910_6_7, eigen_quad_p910_6_8, eigen_quad_p910_6_9, eigen_quad_p910_6_10, eigen_quad_p910_7_2, eigen_quad_p910_7_3, eigen_quad_p910_7_4, eigen_quad_p910_7_5, eigen_quad_p910_7_6, eigen_quad_p910_7_7, eigen_quad_p910_7_8, eigen_quad_p910_7_9, eigen_quad_p910_7_10, eigen_quad_p910_8_2, eigen_quad_p910_8_3, eigen_quad_p910_8_4, eigen_quad_p910_8_5, eigen_quad_p910_8_6, eigen_quad_p910_8_7, eigen_quad_p910_8_8, eigen_quad_p910_8_9, eigen_quad_p910_8_10, eigen_quad_p910_9_2, eigen_quad_p910_9_3, eigen_quad_p910_9_4, eigen_quad_p910_9_5, eigen_quad_p910_9_6, eigen_quad_p910_9_7, eigen_quad_p910_9_8, eigen_quad_p910_9_9, eigen_quad_p910_9_10,

  eigen_quad_p102_2_2, eigen_quad_p102_3_2, eigen_quad_p102_4_2, eigen_quad_p102_5_2, eigen_quad_p102_6_2, eigen_quad_p102_7_2, eigen_quad_p102_8_2, eigen_quad_p102_9_2, eigen_quad_p102_10_2, eigen_quad_p103_2_2, eigen_quad_p103_2_3, eigen_quad_p103_3_2, eigen_quad_p103_3_3, eigen_quad_p103_4_2, eigen_quad_p103_4_3, eigen_quad_p103_5_2, eigen_quad_p103_5_3, eigen_quad_p103_6_2, eigen_quad_p103_6_3, eigen_quad_p103_7_2, eigen_quad_p103_7_3, eigen_quad_p103_8_2, eigen_quad_p103_8_3, eigen_quad_p103_9_2, eigen_quad_p103_9_3, eigen_quad_p103_10_2, eigen_quad_p103_10_3, eigen_quad_p104_2_2, eigen_quad_p104_2_3, eigen_quad_p104_2_4, eigen_quad_p104_3_2, eigen_quad_p104_3_3, eigen_quad_p104_3_4, eigen_quad_p104_4_2, eigen_quad_p104_4_3, eigen_quad_p104_4_4, eigen_quad_p104_5_2, eigen_quad_p104_5_3, eigen_quad_p104_5_4, eigen_quad_p104_6_2, eigen_quad_p104_6_3, eigen_quad_p104_6_4, eigen_quad_p104_7_2, eigen_quad_p104_7_3, eigen_quad_p104_7_4, eigen_quad_p104_8_2, eigen_quad_p104_8_3, eigen_quad_p104_8_4, eigen_quad_p104_9_2, eigen_quad_p104_9_3, eigen_quad_p104_9_4, eigen_quad_p104_10_2, eigen_quad_p104_10_3, eigen_quad_p104_10_4, eigen_quad_p105_2_2, eigen_quad_p105_2_3, eigen_quad_p105_2_4, eigen_quad_p105_2_5, eigen_quad_p105_3_2, eigen_quad_p105_3_3, eigen_quad_p105_3_4, eigen_quad_p105_3_5, eigen_quad_p105_4_2, eigen_quad_p105_4_3, eigen_quad_p105_4_4, eigen_quad_p105_4_5, eigen_quad_p105_5_2, eigen_quad_p105_5_3, eigen_quad_p105_5_4, eigen_quad_p105_5_5, eigen_quad_p105_6_2, eigen_quad_p105_6_3, eigen_quad_p105_6_4, eigen_quad_p105_6_5, eigen_quad_p105_7_2, eigen_quad_p105_7_3, eigen_quad_p105_7_4, eigen_quad_p105_7_5, eigen_quad_p105_8_2, eigen_quad_p105_8_3, eigen_quad_p105_8_4, eigen_quad_p105_8_5, eigen_quad_p105_9_2, eigen_quad_p105_9_3, eigen_quad_p105_9_4, eigen_quad_p105_9_5, eigen_quad_p105_10_2, eigen_quad_p105_10_3, eigen_quad_p105_10_4, eigen_quad_p105_10_5, eigen_quad_p106_2_2, eigen_quad_p106_2_3, eigen_quad_p106_2_4, eigen_quad_p106_2_5, eigen_quad_p106_2_6, eigen_quad_p106_3_2, eigen_quad_p106_3_3, eigen_quad_p106_3_4, eigen_quad_p106_3_5, eigen_quad_p106_3_6, eigen_quad_p106_4_2, eigen_quad_p106_4_3, eigen_quad_p106_4_4, eigen_quad_p106_4_5, eigen_quad_p106_4_6, eigen_quad_p106_5_2, eigen_quad_p106_5_3, eigen_quad_p106_5_4, eigen_quad_p106_5_5, eigen_quad_p106_5_6, eigen_quad_p106_6_2, eigen_quad_p106_6_3, eigen_quad_p106_6_4, eigen_quad_p106_6_5, eigen_quad_p106_6_6, eigen_quad_p106_7_2, eigen_quad_p106_7_3, eigen_quad_p106_7_4, eigen_quad_p106_7_5, eigen_quad_p106_7_6, eigen_quad_p106_8_2, eigen_quad_p106_8_3, eigen_quad_p106_8_4, eigen_quad_p106_8_5, eigen_quad_p106_8_6, eigen_quad_p106_9_2, eigen_quad_p106_9_3, eigen_quad_p106_9_4, eigen_quad_p106_9_5, eigen_quad_p106_9_6, eigen_quad_p106_10_2, eigen_quad_p106_10_3, eigen_quad_p106_10_4, eigen_quad_p106_10_5, eigen_quad_p106_10_6, eigen_quad_p107_2_2, eigen_quad_p107_2_3, eigen_quad_p107_2_4, eigen_quad_p107_2_5, eigen_quad_p107_2_6, eigen_quad_p107_2_7, eigen_quad_p107_3_2, eigen_quad_p107_3_3, eigen_quad_p107_3_4, eigen_quad_p107_3_5, eigen_quad_p107_3_6, eigen_quad_p107_3_7, eigen_quad_p107_4_2, eigen_quad_p107_4_3, eigen_quad_p107_4_4, eigen_quad_p107_4_5, eigen_quad_p107_4_6, eigen_quad_p107_4_7, eigen_quad_p107_5_2, eigen_quad_p107_5_3, eigen_quad_p107_5_4, eigen_quad_p107_5_5, eigen_quad_p107_5_6, eigen_quad_p107_5_7, eigen_quad_p107_6_2, eigen_quad_p107_6_3, eigen_quad_p107_6_4, eigen_quad_p107_6_5, eigen_quad_p107_6_6, eigen_quad_p107_6_7, eigen_quad_p107_7_2, eigen_quad_p107_7_3, eigen_quad_p107_7_4, eigen_quad_p107_7_5, eigen_quad_p107_7_6, eigen_quad_p107_7_7, eigen_quad_p107_8_2, eigen_quad_p107_8_3, eigen_quad_p107_8_4, eigen_quad_p107_8_5, eigen_quad_p107_8_6, eigen_quad_p107_8_7, eigen_quad_p107_9_2, eigen_quad_p107_9_3, eigen_quad_p107_9_4, eigen_quad_p107_9_5, eigen_quad_p107_9_6, eigen_quad_p107_9_7, eigen_quad_p107_10_2, eigen_quad_p107_10_3, eigen_quad_p107_10_4, eigen_quad_p107_10_5, eigen_quad_p107_10_6, eigen_quad_p107_10_7, eigen_quad_p108_2_2, eigen_quad_p108_2_3, eigen_quad_p108_2_4, eigen_quad_p108_2_5, eigen_quad_p108_2_6, eigen_quad_p108_2_7, eigen_quad_p108_2_8, eigen_quad_p108_3_2, eigen_quad_p108_3_3, eigen_quad_p108_3_4, eigen_quad_p108_3_5, eigen_quad_p108_3_6, eigen_quad_p108_3_7, eigen_quad_p108_3_8, eigen_quad_p108_4_2, eigen_quad_p108_4_3, eigen_quad_p108_4_4, eigen_quad_p108_4_5, eigen_quad_p108_4_6, eigen_quad_p108_4_7, eigen_quad_p108_4_8, eigen_quad_p108_5_2, eigen_quad_p108_5_3, eigen_quad_p108_5_4, eigen_quad_p108_5_5, eigen_quad_p108_5_6, eigen_quad_p108_5_7, eigen_quad_p108_5_8, eigen_quad_p108_6_2, eigen_quad_p108_6_3, eigen_quad_p108_6_4, eigen_quad_p108_6_5, eigen_quad_p108_6_6, eigen_quad_p108_6_7, eigen_quad_p108_6_8, eigen_quad_p108_7_2, eigen_quad_p108_7_3, eigen_quad_p108_7_4, eigen_quad_p108_7_5, eigen_quad_p108_7_6, eigen_quad_p108_7_7, eigen_quad_p108_7_8, eigen_quad_p108_8_2, eigen_quad_p108_8_3, eigen_quad_p108_8_4, eigen_quad_p108_8_5, eigen_quad_p108_8_6, eigen_quad_p108_8_7, eigen_quad_p108_8_8, eigen_quad_p108_9_2, eigen_quad_p108_9_3, eigen_quad_p108_9_4, eigen_quad_p108_9_5, eigen_quad_p108_9_6, eigen_quad_p108_9_7, eigen_quad_p108_9_8, eigen_quad_p108_10_2, eigen_quad_p108_10_3, eigen_quad_p108_10_4, eigen_quad_p108_10_5, eigen_quad_p108_10_6, eigen_quad_p108_10_7, eigen_quad_p108_10_8, eigen_quad_p109_2_2, eigen_quad_p109_2_3, eigen_quad_p109_2_4, eigen_quad_p109_2_5, eigen_quad_p109_2_6, eigen_quad_p109_2_7, eigen_quad_p109_2_8, eigen_quad_p109_2_9, eigen_quad_p109_3_2, eigen_quad_p109_3_3, eigen_quad_p109_3_4, eigen_quad_p109_3_5, eigen_quad_p109_3_6, eigen_quad_p109_3_7, eigen_quad_p109_3_8, eigen_quad_p109_3_9, eigen_quad_p109_4_2, eigen_quad_p109_4_3, eigen_quad_p109_4_4, eigen_quad_p109_4_5, eigen_quad_p109_4_6, eigen_quad_p109_4_7, eigen_quad_p109_4_8, eigen_quad_p109_4_9, eigen_quad_p109_5_2, eigen_quad_p109_5_3, eigen_quad_p109_5_4, eigen_quad_p109_5_5, eigen_quad_p109_5_6, eigen_quad_p109_5_7, eigen_quad_p109_5_8, eigen_quad_p109_5_9, eigen_quad_p109_6_2, eigen_quad_p109_6_3, eigen_quad_p109_6_4, eigen_quad_p109_6_5, eigen_quad_p109_6_6, eigen_quad_p109_6_7, eigen_quad_p109_6_8, eigen_quad_p109_6_9, eigen_quad_p109_7_2, eigen_quad_p109_7_3, eigen_quad_p109_7_4, eigen_quad_p109_7_5, eigen_quad_p109_7_6, eigen_quad_p109_7_7, eigen_quad_p109_7_8, eigen_quad_p109_7_9, eigen_quad_p109_8_2, eigen_quad_p109_8_3, eigen_quad_p109_8_4, eigen_quad_p109_8_5, eigen_quad_p109_8_6, eigen_quad_p109_8_7, eigen_quad_p109_8_8, eigen_quad_p109_8_9, eigen_quad_p109_9_2, eigen_quad_p109_9_3, eigen_quad_p109_9_4, eigen_quad_p109_9_5, eigen_quad_p109_9_6, eigen_quad_p109_9_7, eigen_quad_p109_9_8, eigen_quad_p109_9_9, eigen_quad_p109_10_2, eigen_quad_p109_10_3, eigen_quad_p109_10_4, eigen_quad_p109_10_5, eigen_quad_p109_10_6, eigen_quad_p109_10_7, eigen_quad_p109_10_8, eigen_quad_p109_10_9, eigen_quad_p1010_2_2, eigen_quad_p1010_2_3, eigen_quad_p1010_2_4, eigen_quad_p1010_2_5, eigen_quad_p1010_2_6, eigen_quad_p1010_2_7, eigen_quad_p1010_2_8, eigen_quad_p1010_2_9, eigen_quad_p1010_2_10, eigen_quad_p1010_3_2, eigen_quad_p1010_3_3, eigen_quad_p1010_3_4, eigen_quad_p1010_3_5, eigen_quad_p1010_3_6, eigen_quad_p1010_3_7, eigen_quad_p1010_3_8, eigen_quad_p1010_3_9, eigen_quad_p1010_3_10, eigen_quad_p1010_4_2, eigen_quad_p1010_4_3, eigen_quad_p1010_4_4, eigen_quad_p1010_4_5, eigen_quad_p1010_4_6, eigen_quad_p1010_4_7, eigen_quad_p1010_4_8, eigen_quad_p1010_4_9, eigen_quad_p1010_4_10, eigen_quad_p1010_5_2, eigen_quad_p1010_5_3, eigen_quad_p1010_5_4, eigen_quad_p1010_5_5, eigen_quad_p1010_5_6, eigen_quad_p1010_5_7, eigen_quad_p1010_5_8, eigen_quad_p1010_5_9, eigen_quad_p1010_5_10, eigen_quad_p1010_6_2, eigen_quad_p1010_6_3, eigen_quad_p1010_6_4, eigen_quad_p1010_6_5, eigen_quad_p1010_6_6, eigen_quad_p1010_6_7, eigen_quad_p1010_6_8, eigen_quad_p1010_6_9, eigen_quad_p1010_6_10, eigen_quad_p1010_7_2, eigen_quad_p1010_7_3, eigen_quad_p1010_7_4, eigen_quad_p1010_7_5, eigen_quad_p1010_7_6, eigen_quad_p1010_7_7, eigen_quad_p1010_7_8, eigen_quad_p1010_7_9, eigen_quad_p1010_7_10, eigen_quad_p1010_8_2, eigen_quad_p1010_8_3, eigen_quad_p1010_8_4, eigen_quad_p1010_8_5, eigen_quad_p1010_8_6, eigen_quad_p1010_8_7, eigen_quad_p1010_8_8, eigen_quad_p1010_8_9, eigen_quad_p1010_8_10, eigen_quad_p1010_9_2, eigen_quad_p1010_9_3, eigen_quad_p1010_9_4, eigen_quad_p1010_9_5, eigen_quad_p1010_9_6, eigen_quad_p1010_9_7, eigen_quad_p1010_9_8, eigen_quad_p1010_9_9, eigen_quad_p1010_9_10, eigen_quad_p1010_10_2, eigen_quad_p1010_10_3, eigen_quad_p1010_10_4, eigen_quad_p1010_10_5, eigen_quad_p1010_10_6, eigen_quad_p1010_10_7, eigen_quad_p1010_10_8, eigen_quad_p1010_10_9, eigen_quad_p1010_10_10
};


static Shapeset::shape_fn_t eigen_quad_fn_dx[] =
{
  eigen_quad_0_0x,   eigen_quad_0_1x,   eigen_quad_0_2x,   eigen_quad_0_3x_0, eigen_quad_0_3x_1,  eigen_quad_0_4x,   eigen_quad_0_5x_0,  eigen_quad_0_5x_1, eigen_quad_0_6x,   eigen_quad_0_7x_0, eigen_quad_0_7x_1, eigen_quad_0_8x,  eigen_quad_0_9x_0,  eigen_quad_0_9x_1, eigen_quad_0_10x,  eigen_quad_1_0x,  eigen_quad_1_1x,   eigen_quad_1_2x,   eigen_quad_1_3x_0, eigen_quad_1_3x_1, eigen_quad_1_4x,    eigen_quad_1_5x_0, eigen_quad_1_5x_1, eigen_quad_1_6x,  eigen_quad_1_7x_0, eigen_quad_1_7x_1, eigen_quad_1_8x,   eigen_quad_1_9x_0, eigen_quad_1_9x_1,  eigen_quad_1_10x,  eigen_quad_2_0x,   eigen_quad_2_1x,   eigen_quad_3_0x_0, eigen_quad_3_0x_1, eigen_quad_3_1x_0, eigen_quad_3_1x_1, eigen_quad_4_0x,   eigen_quad_4_1x,   eigen_quad_5_0x_0, eigen_quad_5_0x_1, eigen_quad_5_1x_0, eigen_quad_5_1x_1, eigen_quad_6_0x,   eigen_quad_6_1x,   eigen_quad_7_0x_0, eigen_quad_7_0x_1, eigen_quad_7_1x_0, eigen_quad_7_1x_1, eigen_quad_8_0x,   eigen_quad_8_1x,   eigen_quad_9_0x_0, eigen_quad_9_0x_1, eigen_quad_9_1x_0, eigen_quad_9_1x_1, eigen_quad_10_0x,  eigen_quad_10_1x,

  eigen_quad_p22_2_2x, eigen_quad_p23_2_2x, eigen_quad_p23_2_3x, eigen_quad_p24_2_2x, eigen_quad_p24_2_3x, eigen_quad_p24_2_4x, eigen_quad_p25_2_2x, eigen_quad_p25_2_3x, eigen_quad_p25_2_4x, eigen_quad_p25_2_5x, eigen_quad_p26_2_2x, eigen_quad_p26_2_3x, eigen_quad_p26_2_4x, eigen_quad_p26_2_5x, eigen_quad_p26_2_6x, eigen_quad_p27_2_2x, eigen_quad_p27_2_3x, eigen_quad_p27_2_4x, eigen_quad_p27_2_5x, eigen_quad_p27_2_6x, eigen_quad_p27_2_7x, eigen_quad_p28_2_2x, eigen_quad_p28_2_3x, eigen_quad_p28_2_4x, eigen_quad_p28_2_5x, eigen_quad_p28_2_6x, eigen_quad_p28_2_7x, eigen_quad_p28_2_8x, eigen_quad_p29_2_2x, eigen_quad_p29_2_3x, eigen_quad_p29_2_4x, eigen_quad_p29_2_5x, eigen_quad_p29_2_6x, eigen_quad_p29_2_7x, eigen_quad_p29_2_8x, eigen_quad_p29_2_9x, eigen_quad_p210_2_2x, eigen_quad_p210_2_3x, eigen_quad_p210_2_4x, eigen_quad_p210_2_5x, eigen_quad_p210_2_6x, eigen_quad_p210_2_7x, eigen_quad_p210_2_8x, eigen_quad_p210_2_9x, eigen_quad_p210_2_10x,

  eigen_quad_p32_2_2x, eigen_quad_p32_3_2x, eigen_quad_p33_2_2x, eigen_quad_p33_2_3x, eigen_quad_p33_3_2x, eigen_quad_p33_3_3x, eigen_quad_p34_2_2x, eigen_quad_p34_2_3x, eigen_quad_p34_2_4x, eigen_quad_p34_3_2x, eigen_quad_p34_3_3x, eigen_quad_p34_3_4x, eigen_quad_p35_2_2x, eigen_quad_p35_2_3x, eigen_quad_p35_2_4x, eigen_quad_p35_2_5x, eigen_quad_p35_3_2x, eigen_quad_p35_3_3x, eigen_quad_p35_3_4x, eigen_quad_p35_3_5x, eigen_quad_p36_2_2x, eigen_quad_p36_2_3x, eigen_quad_p36_2_4x, eigen_quad_p36_2_5x, eigen_quad_p36_2_6x, eigen_quad_p36_3_2x, eigen_quad_p36_3_3x, eigen_quad_p36_3_4x, eigen_quad_p36_3_5x, eigen_quad_p36_3_6x, eigen_quad_p37_2_2x, eigen_quad_p37_2_3x, eigen_quad_p37_2_4x, eigen_quad_p37_2_5x, eigen_quad_p37_2_6x, eigen_quad_p37_2_7x, eigen_quad_p37_3_2x, eigen_quad_p37_3_3x, eigen_quad_p37_3_4x, eigen_quad_p37_3_5x, eigen_quad_p37_3_6x, eigen_quad_p37_3_7x, eigen_quad_p38_2_2x, eigen_quad_p38_2_3x, eigen_quad_p38_2_4x, eigen_quad_p38_2_5x, eigen_quad_p38_2_6x, eigen_quad_p38_2_7x, eigen_quad_p38_2_8x, eigen_quad_p38_3_2x, eigen_quad_p38_3_3x, eigen_quad_p38_3_4x, eigen_quad_p38_3_5x, eigen_quad_p38_3_6x, eigen_quad_p38_3_7x, eigen_quad_p38_3_8x, eigen_quad_p39_2_2x, eigen_quad_p39_2_3x, eigen_quad_p39_2_4x, eigen_quad_p39_2_5x, eigen_quad_p39_2_6x, eigen_quad_p39_2_7x, eigen_quad_p39_2_8x, eigen_quad_p39_2_9x, eigen_quad_p39_3_2x, eigen_quad_p39_3_3x, eigen_quad_p39_3_4x, eigen_quad_p39_3_5x, eigen_quad_p39_3_6x, eigen_quad_p39_3_7x, eigen_quad_p39_3_8x, eigen_quad_p39_3_9x, eigen_quad_p310_2_2x, eigen_quad_p310_2_3x, eigen_quad_p310_2_4x, eigen_quad_p310_2_5x, eigen_quad_p310_2_6x, eigen_quad_p310_2_7x, eigen_quad_p310_2_8x, eigen_quad_p310_2_9x, eigen_quad_p310_2_10x, eigen_quad_p310_3_2x, eigen_quad_p310_3_3x, eigen_quad_p310_3_4x, eigen_quad_p310_3_5x, eigen_quad_p310_3_6x, eigen_quad_p310_3_7x, eigen_quad_p310_3_8x, eigen_quad_p310_3_9x, eigen_quad_p310_3_10x,

  eigen_quad_p42_2_2x, eigen_quad_p42_3_2x, eigen_quad_p42_4_2x, eigen_quad_p43_2_2x, eigen_quad_p43_2_3x, eigen_quad_p43_3_2x, eigen_quad_p43_3_3x, eigen_quad_p43_4_2x, eigen_quad_p43_4_3x, eigen_quad_p44_2_2x, eigen_quad_p44_2_3x, eigen_quad_p44_2_4x, eigen_quad_p44_3_2x, eigen_quad_p44_3_3x, eigen_quad_p44_3_4x, eigen_quad_p44_4_2x, eigen_quad_p44_4_3x, eigen_quad_p44_4_4x, eigen_quad_p45_2_2x, eigen_quad_p45_2_3x, eigen_quad_p45_2_4x, eigen_quad_p45_2_5x, eigen_quad_p45_3_2x, eigen_quad_p45_3_3x, eigen_quad_p45_3_4x, eigen_quad_p45_3_5x, eigen_quad_p45_4_2x, eigen_quad_p45_4_3x, eigen_quad_p45_4_4x, eigen_quad_p45_4_5x, eigen_quad_p46_2_2x, eigen_quad_p46_2_3x, eigen_quad_p46_2_4x, eigen_quad_p46_2_5x, eigen_quad_p46_2_6x, eigen_quad_p46_3_2x, eigen_quad_p46_3_3x, eigen_quad_p46_3_4x, eigen_quad_p46_3_5x, eigen_quad_p46_3_6x, eigen_quad_p46_4_2x, eigen_quad_p46_4_3x, eigen_quad_p46_4_4x, eigen_quad_p46_4_5x, eigen_quad_p46_4_6x, eigen_quad_p47_2_2x, eigen_quad_p47_2_3x, eigen_quad_p47_2_4x, eigen_quad_p47_2_5x, eigen_quad_p47_2_6x, eigen_quad_p47_2_7x, eigen_quad_p47_3_2x, eigen_quad_p47_3_3x, eigen_quad_p47_3_4x, eigen_quad_p47_3_5x, eigen_quad_p47_3_6x, eigen_quad_p47_3_7x, eigen_quad_p47_4_2x, eigen_quad_p47_4_3x, eigen_quad_p47_4_4x, eigen_quad_p47_4_5x, eigen_quad_p47_4_6x, eigen_quad_p47_4_7x, eigen_quad_p48_2_2x, eigen_quad_p48_2_3x, eigen_quad_p48_2_4x, eigen_quad_p48_2_5x, eigen_quad_p48_2_6x, eigen_quad_p48_2_7x, eigen_quad_p48_2_8x, eigen_quad_p48_3_2x, eigen_quad_p48_3_3x, eigen_quad_p48_3_4x, eigen_quad_p48_3_5x, eigen_quad_p48_3_6x, eigen_quad_p48_3_7x, eigen_quad_p48_3_8x, eigen_quad_p48_4_2x, eigen_quad_p48_4_3x, eigen_quad_p48_4_4x, eigen_quad_p48_4_5x, eigen_quad_p48_4_6x, eigen_quad_p48_4_7x, eigen_quad_p48_4_8x, eigen_quad_p49_2_2x, eigen_quad_p49_2_3x, eigen_quad_p49_2_4x, eigen_quad_p49_2_5x, eigen_quad_p49_2_6x, eigen_quad_p49_2_7x, eigen_quad_p49_2_8x, eigen_quad_p49_2_9x, eigen_quad_p49_3_2x, eigen_quad_p49_3_3x, eigen_quad_p49_3_4x, eigen_quad_p49_3_5x, eigen_quad_p49_3_6x, eigen_quad_p49_3_7x, eigen_quad_p49_3_8x, eigen_quad_p49_3_9x, eigen_quad_p49_4_2x, eigen_quad_p49_4_3x, eigen_quad_p49_4_4x, eigen_quad_p49_4_5x, eigen_quad_p49_4_6x, eigen_quad_p49_4_7x, eigen_quad_p49_4_8x, eigen_quad_p49_4_9x, eigen_quad_p410_2_2x, eigen_quad_p410_2_3x, eigen_quad_p410_2_4x, eigen_quad_p410_2_5x, eigen_quad_p410_2_6x, eigen_quad_p410_2_7x, eigen_quad_p410_2_8x, eigen_quad_p410_2_9x, eigen_quad_p410_2_10x, eigen_quad_p410_3_2x, eigen_quad_p410_3_3x, eigen_quad_p410_3_4x, eigen_quad_p410_3_5x, eigen_quad_p410_3_6x, eigen_quad_p410_3_7x, eigen_quad_p410_3_8x, eigen_quad_p410_3_9x, eigen_quad_p410_3_10x, eigen_quad_p410_4_2x, eigen_quad_p410_4_3x, eigen_quad_p410_4_4x, eigen_quad_p410_4_5x, eigen_quad_p410_4_6x, eigen_quad_p410_4_7x, eigen_quad_p410_4_8x, eigen_quad_p410_4_9x, eigen_quad_p410_4_10x,

  eigen_quad_p52_2_2x, eigen_quad_p52_3_2x, eigen_quad_p52_4_2x, eigen_quad_p52_5_2x, eigen_quad_p53_2_2x, eigen_quad_p53_2_3x, eigen_quad_p53_3_2x, eigen_quad_p53_3_3x, eigen_quad_p53_4_2x, eigen_quad_p53_4_3x, eigen_quad_p53_5_2x, eigen_quad_p53_5_3x, eigen_quad_p54_2_2x, eigen_quad_p54_2_3x, eigen_quad_p54_2_4x, eigen_quad_p54_3_2x, eigen_quad_p54_3_3x, eigen_quad_p54_3_4x, eigen_quad_p54_4_2x, eigen_quad_p54_4_3x, eigen_quad_p54_4_4x, eigen_quad_p54_5_2x, eigen_quad_p54_5_3x, eigen_quad_p54_5_4x, eigen_quad_p55_2_2x, eigen_quad_p55_2_3x, eigen_quad_p55_2_4x, eigen_quad_p55_2_5x, eigen_quad_p55_3_2x, eigen_quad_p55_3_3x, eigen_quad_p55_3_4x, eigen_quad_p55_3_5x, eigen_quad_p55_4_2x, eigen_quad_p55_4_3x, eigen_quad_p55_4_4x, eigen_quad_p55_4_5x, eigen_quad_p55_5_2x, eigen_quad_p55_5_3x, eigen_quad_p55_5_4x, eigen_quad_p55_5_5x, eigen_quad_p56_2_2x, eigen_quad_p56_2_3x, eigen_quad_p56_2_4x, eigen_quad_p56_2_5x, eigen_quad_p56_2_6x, eigen_quad_p56_3_2x, eigen_quad_p56_3_3x, eigen_quad_p56_3_4x, eigen_quad_p56_3_5x, eigen_quad_p56_3_6x, eigen_quad_p56_4_2x, eigen_quad_p56_4_3x, eigen_quad_p56_4_4x, eigen_quad_p56_4_5x, eigen_quad_p56_4_6x, eigen_quad_p56_5_2x, eigen_quad_p56_5_3x, eigen_quad_p56_5_4x, eigen_quad_p56_5_5x, eigen_quad_p56_5_6x, eigen_quad_p57_2_2x, eigen_quad_p57_2_3x, eigen_quad_p57_2_4x, eigen_quad_p57_2_5x, eigen_quad_p57_2_6x, eigen_quad_p57_2_7x, eigen_quad_p57_3_2x, eigen_quad_p57_3_3x, eigen_quad_p57_3_4x, eigen_quad_p57_3_5x, eigen_quad_p57_3_6x, eigen_quad_p57_3_7x, eigen_quad_p57_4_2x, eigen_quad_p57_4_3x, eigen_quad_p57_4_4x, eigen_quad_p57_4_5x, eigen_quad_p57_4_6x, eigen_quad_p57_4_7x, eigen_quad_p57_5_2x, eigen_quad_p57_5_3x, eigen_quad_p57_5_4x, eigen_quad_p57_5_5x, eigen_quad_p57_5_6x, eigen_quad_p57_5_7x, eigen_quad_p58_2_2x, eigen_quad_p58_2_3x, eigen_quad_p58_2_4x, eigen_quad_p58_2_5x, eigen_quad_p58_2_6x, eigen_quad_p58_2_7x, eigen_quad_p58_2_8x, eigen_quad_p58_3_2x, eigen_quad_p58_3_3x, eigen_quad_p58_3_4x, eigen_quad_p58_3_5x, eigen_quad_p58_3_6x, eigen_quad_p58_3_7x, eigen_quad_p58_3_8x, eigen_quad_p58_4_2x, eigen_quad_p58_4_3x, eigen_quad_p58_4_4x, eigen_quad_p58_4_5x, eigen_quad_p58_4_6x, eigen_quad_p58_4_7x, eigen_quad_p58_4_8x, eigen_quad_p58_5_2x, eigen_quad_p58_5_3x, eigen_quad_p58_5_4x, eigen_quad_p58_5_5x, eigen_quad_p58_5_6x, eigen_quad_p58_5_7x, eigen_quad_p58_5_8x, eigen_quad_p59_2_2x, eigen_quad_p59_2_3x, eigen_quad_p59_2_4x, eigen_quad_p59_2_5x, eigen_quad_p59_2_6x, eigen_quad_p59_2_7x, eigen_quad_p59_2_8x, eigen_quad_p59_2_9x, eigen_quad_p59_3_2x, eigen_quad_p59_3_3x, eigen_quad_p59_3_4x, eigen_quad_p59_3_5x, eigen_quad_p59_3_6x, eigen_quad_p59_3_7x, eigen_quad_p59_3_8x, eigen_quad_p59_3_9x, eigen_quad_p59_4_2x, eigen_quad_p59_4_3x, eigen_quad_p59_4_4x, eigen_quad_p59_4_5x, eigen_quad_p59_4_6x, eigen_quad_p59_4_7x, eigen_quad_p59_4_8x, eigen_quad_p59_4_9x, eigen_quad_p59_5_2x, eigen_quad_p59_5_3x, eigen_quad_p59_5_4x, eigen_quad_p59_5_5x, eigen_quad_p59_5_6x, eigen_quad_p59_5_7x, eigen_quad_p59_5_8x, eigen_quad_p59_5_9x, eigen_quad_p510_2_2x, eigen_quad_p510_2_3x, eigen_quad_p510_2_4x, eigen_quad_p510_2_5x, eigen_quad_p510_2_6x, eigen_quad_p510_2_7x, eigen_quad_p510_2_8x, eigen_quad_p510_2_9x, eigen_quad_p510_2_10x, eigen_quad_p510_3_2x, eigen_quad_p510_3_3x, eigen_quad_p510_3_4x, eigen_quad_p510_3_5x, eigen_quad_p510_3_6x, eigen_quad_p510_3_7x, eigen_quad_p510_3_8x, eigen_quad_p510_3_9x, eigen_quad_p510_3_10x, eigen_quad_p510_4_2x, eigen_quad_p510_4_3x, eigen_quad_p510_4_4x, eigen_quad_p510_4_5x, eigen_quad_p510_4_6x, eigen_quad_p510_4_7x, eigen_quad_p510_4_8x, eigen_quad_p510_4_9x, eigen_quad_p510_4_10x, eigen_quad_p510_5_2x, eigen_quad_p510_5_3x, eigen_quad_p510_5_4x, eigen_quad_p510_5_5x, eigen_quad_p510_5_6x, eigen_quad_p510_5_7x, eigen_quad_p510_5_8x, eigen_quad_p510_5_9x, eigen_quad_p510_5_10x,

  eigen_quad_p62_2_2x, eigen_quad_p62_3_2x, eigen_quad_p62_4_2x, eigen_quad_p62_5_2x, eigen_quad_p62_6_2x, eigen_quad_p63_2_2x, eigen_quad_p63_2_3x, eigen_quad_p63_3_2x, eigen_quad_p63_3_3x, eigen_quad_p63_4_2x, eigen_quad_p63_4_3x, eigen_quad_p63_5_2x, eigen_quad_p63_5_3x, eigen_quad_p63_6_2x, eigen_quad_p63_6_3x, eigen_quad_p64_2_2x, eigen_quad_p64_2_3x, eigen_quad_p64_2_4x, eigen_quad_p64_3_2x, eigen_quad_p64_3_3x, eigen_quad_p64_3_4x, eigen_quad_p64_4_2x, eigen_quad_p64_4_3x, eigen_quad_p64_4_4x, eigen_quad_p64_5_2x, eigen_quad_p64_5_3x, eigen_quad_p64_5_4x, eigen_quad_p64_6_2x, eigen_quad_p64_6_3x, eigen_quad_p64_6_4x, eigen_quad_p65_2_2x, eigen_quad_p65_2_3x, eigen_quad_p65_2_4x, eigen_quad_p65_2_5x, eigen_quad_p65_3_2x, eigen_quad_p65_3_3x, eigen_quad_p65_3_4x, eigen_quad_p65_3_5x, eigen_quad_p65_4_2x, eigen_quad_p65_4_3x, eigen_quad_p65_4_4x, eigen_quad_p65_4_5x, eigen_quad_p65_5_2x, eigen_quad_p65_5_3x, eigen_quad_p65_5_4x, eigen_quad_p65_5_5x, eigen_quad_p65_6_2x, eigen_quad_p65_6_3x, eigen_quad_p65_6_4x, eigen_quad_p65_6_5x, eigen_quad_p66_2_2x, eigen_quad_p66_2_3x, eigen_quad_p66_2_4x, eigen_quad_p66_2_5x, eigen_quad_p66_2_6x, eigen_quad_p66_3_2x, eigen_quad_p66_3_3x, eigen_quad_p66_3_4x, eigen_quad_p66_3_5x, eigen_quad_p66_3_6x, eigen_quad_p66_4_2x, eigen_quad_p66_4_3x, eigen_quad_p66_4_4x, eigen_quad_p66_4_5x, eigen_quad_p66_4_6x, eigen_quad_p66_5_2x, eigen_quad_p66_5_3x, eigen_quad_p66_5_4x, eigen_quad_p66_5_5x, eigen_quad_p66_5_6x, eigen_quad_p66_6_2x, eigen_quad_p66_6_3x, eigen_quad_p66_6_4x, eigen_quad_p66_6_5x, eigen_quad_p66_6_6x, eigen_quad_p67_2_2x, eigen_quad_p67_2_3x, eigen_quad_p67_2_4x, eigen_quad_p67_2_5x, eigen_quad_p67_2_6x, eigen_quad_p67_2_7x, eigen_quad_p67_3_2x, eigen_quad_p67_3_3x, eigen_quad_p67_3_4x, eigen_quad_p67_3_5x, eigen_quad_p67_3_6x, eigen_quad_p67_3_7x, eigen_quad_p67_4_2x, eigen_quad_p67_4_3x, eigen_quad_p67_4_4x, eigen_quad_p67_4_5x, eigen_quad_p67_4_6x, eigen_quad_p67_4_7x, eigen_quad_p67_5_2x, eigen_quad_p67_5_3x, eigen_quad_p67_5_4x, eigen_quad_p67_5_5x, eigen_quad_p67_5_6x, eigen_quad_p67_5_7x, eigen_quad_p67_6_2x, eigen_quad_p67_6_3x, eigen_quad_p67_6_4x, eigen_quad_p67_6_5x, eigen_quad_p67_6_6x, eigen_quad_p67_6_7x, eigen_quad_p68_2_2x, eigen_quad_p68_2_3x, eigen_quad_p68_2_4x, eigen_quad_p68_2_5x, eigen_quad_p68_2_6x, eigen_quad_p68_2_7x, eigen_quad_p68_2_8x, eigen_quad_p68_3_2x, eigen_quad_p68_3_3x, eigen_quad_p68_3_4x, eigen_quad_p68_3_5x, eigen_quad_p68_3_6x, eigen_quad_p68_3_7x, eigen_quad_p68_3_8x, eigen_quad_p68_4_2x, eigen_quad_p68_4_3x, eigen_quad_p68_4_4x, eigen_quad_p68_4_5x, eigen_quad_p68_4_6x, eigen_quad_p68_4_7x, eigen_quad_p68_4_8x, eigen_quad_p68_5_2x, eigen_quad_p68_5_3x, eigen_quad_p68_5_4x, eigen_quad_p68_5_5x, eigen_quad_p68_5_6x, eigen_quad_p68_5_7x, eigen_quad_p68_5_8x, eigen_quad_p68_6_2x, eigen_quad_p68_6_3x, eigen_quad_p68_6_4x, eigen_quad_p68_6_5x, eigen_quad_p68_6_6x, eigen_quad_p68_6_7x, eigen_quad_p68_6_8x, eigen_quad_p69_2_2x, eigen_quad_p69_2_3x, eigen_quad_p69_2_4x, eigen_quad_p69_2_5x, eigen_quad_p69_2_6x, eigen_quad_p69_2_7x, eigen_quad_p69_2_8x, eigen_quad_p69_2_9x, eigen_quad_p69_3_2x, eigen_quad_p69_3_3x, eigen_quad_p69_3_4x, eigen_quad_p69_3_5x, eigen_quad_p69_3_6x, eigen_quad_p69_3_7x, eigen_quad_p69_3_8x, eigen_quad_p69_3_9x, eigen_quad_p69_4_2x, eigen_quad_p69_4_3x, eigen_quad_p69_4_4x, eigen_quad_p69_4_5x, eigen_quad_p69_4_6x, eigen_quad_p69_4_7x, eigen_quad_p69_4_8x, eigen_quad_p69_4_9x, eigen_quad_p69_5_2x, eigen_quad_p69_5_3x, eigen_quad_p69_5_4x, eigen_quad_p69_5_5x, eigen_quad_p69_5_6x, eigen_quad_p69_5_7x, eigen_quad_p69_5_8x, eigen_quad_p69_5_9x, eigen_quad_p69_6_2x, eigen_quad_p69_6_3x, eigen_quad_p69_6_4x, eigen_quad_p69_6_5x, eigen_quad_p69_6_6x, eigen_quad_p69_6_7x, eigen_quad_p69_6_8x, eigen_quad_p69_6_9x, eigen_quad_p610_2_2x, eigen_quad_p610_2_3x, eigen_quad_p610_2_4x, eigen_quad_p610_2_5x, eigen_quad_p610_2_6x, eigen_quad_p610_2_7x, eigen_quad_p610_2_8x, eigen_quad_p610_2_9x, eigen_quad_p610_2_10x, eigen_quad_p610_3_2x, eigen_quad_p610_3_3x, eigen_quad_p610_3_4x, eigen_quad_p610_3_5x, eigen_quad_p610_3_6x, eigen_quad_p610_3_7x, eigen_quad_p610_3_8x, eigen_quad_p610_3_9x, eigen_quad_p610_3_10x, eigen_quad_p610_4_2x, eigen_quad_p610_4_3x, eigen_quad_p610_4_4x, eigen_quad_p610_4_5x, eigen_quad_p610_4_6x, eigen_quad_p610_4_7x, eigen_quad_p610_4_8x, eigen_quad_p610_4_9x, eigen_quad_p610_4_10x, eigen_quad_p610_5_2x, eigen_quad_p610_5_3x, eigen_quad_p610_5_4x, eigen_quad_p610_5_5x, eigen_quad_p610_5_6x, eigen_quad_p610_5_7x, eigen_quad_p610_5_8x, eigen_quad_p610_5_9x, eigen_quad_p610_5_10x, eigen_quad_p610_6_2x, eigen_quad_p610_6_3x, eigen_quad_p610_6_4x, eigen_quad_p610_6_5x, eigen_quad_p610_6_6x, eigen_quad_p610_6_7x, eigen_quad_p610_6_8x, eigen_quad_p610_6_9x, eigen_quad_p610_6_10x,

  eigen_quad_p72_2_2x, eigen_quad_p72_3_2x, eigen_quad_p72_4_2x, eigen_quad_p72_5_2x, eigen_quad_p72_6_2x, eigen_quad_p72_7_2x, eigen_quad_p73_2_2x, eigen_quad_p73_2_3x, eigen_quad_p73_3_2x, eigen_quad_p73_3_3x, eigen_quad_p73_4_2x, eigen_quad_p73_4_3x, eigen_quad_p73_5_2x, eigen_quad_p73_5_3x, eigen_quad_p73_6_2x, eigen_quad_p73_6_3x, eigen_quad_p73_7_2x, eigen_quad_p73_7_3x, eigen_quad_p74_2_2x, eigen_quad_p74_2_3x, eigen_quad_p74_2_4x, eigen_quad_p74_3_2x, eigen_quad_p74_3_3x, eigen_quad_p74_3_4x, eigen_quad_p74_4_2x, eigen_quad_p74_4_3x, eigen_quad_p74_4_4x, eigen_quad_p74_5_2x, eigen_quad_p74_5_3x, eigen_quad_p74_5_4x, eigen_quad_p74_6_2x, eigen_quad_p74_6_3x, eigen_quad_p74_6_4x, eigen_quad_p74_7_2x, eigen_quad_p74_7_3x, eigen_quad_p74_7_4x, eigen_quad_p75_2_2x, eigen_quad_p75_2_3x, eigen_quad_p75_2_4x, eigen_quad_p75_2_5x, eigen_quad_p75_3_2x, eigen_quad_p75_3_3x, eigen_quad_p75_3_4x, eigen_quad_p75_3_5x, eigen_quad_p75_4_2x, eigen_quad_p75_4_3x, eigen_quad_p75_4_4x, eigen_quad_p75_4_5x, eigen_quad_p75_5_2x, eigen_quad_p75_5_3x, eigen_quad_p75_5_4x, eigen_quad_p75_5_5x, eigen_quad_p75_6_2x, eigen_quad_p75_6_3x, eigen_quad_p75_6_4x, eigen_quad_p75_6_5x, eigen_quad_p75_7_2x, eigen_quad_p75_7_3x, eigen_quad_p75_7_4x, eigen_quad_p75_7_5x, eigen_quad_p76_2_2x, eigen_quad_p76_2_3x, eigen_quad_p76_2_4x, eigen_quad_p76_2_5x, eigen_quad_p76_2_6x, eigen_quad_p76_3_2x, eigen_quad_p76_3_3x, eigen_quad_p76_3_4x, eigen_quad_p76_3_5x, eigen_quad_p76_3_6x, eigen_quad_p76_4_2x, eigen_quad_p76_4_3x, eigen_quad_p76_4_4x, eigen_quad_p76_4_5x, eigen_quad_p76_4_6x, eigen_quad_p76_5_2x, eigen_quad_p76_5_3x, eigen_quad_p76_5_4x, eigen_quad_p76_5_5x, eigen_quad_p76_5_6x, eigen_quad_p76_6_2x, eigen_quad_p76_6_3x, eigen_quad_p76_6_4x, eigen_quad_p76_6_5x, eigen_quad_p76_6_6x, eigen_quad_p76_7_2x, eigen_quad_p76_7_3x, eigen_quad_p76_7_4x, eigen_quad_p76_7_5x, eigen_quad_p76_7_6x, eigen_quad_p77_2_2x, eigen_quad_p77_2_3x, eigen_quad_p77_2_4x, eigen_quad_p77_2_5x, eigen_quad_p77_2_6x, eigen_quad_p77_2_7x, eigen_quad_p77_3_2x, eigen_quad_p77_3_3x, eigen_quad_p77_3_4x, eigen_quad_p77_3_5x, eigen_quad_p77_3_6x, eigen_quad_p77_3_7x, eigen_quad_p77_4_2x, eigen_quad_p77_4_3x, eigen_quad_p77_4_4x, eigen_quad_p77_4_5x, eigen_quad_p77_4_6x, eigen_quad_p77_4_7x, eigen_quad_p77_5_2x, eigen_quad_p77_5_3x, eigen_quad_p77_5_4x, eigen_quad_p77_5_5x, eigen_quad_p77_5_6x, eigen_quad_p77_5_7x, eigen_quad_p77_6_2x, eigen_quad_p77_6_3x, eigen_quad_p77_6_4x, eigen_quad_p77_6_5x, eigen_quad_p77_6_6x, eigen_quad_p77_6_7x, eigen_quad_p77_7_2x, eigen_quad_p77_7_3x, eigen_quad_p77_7_4x, eigen_quad_p77_7_5x, eigen_quad_p77_7_6x, eigen_quad_p77_7_7x, eigen_quad_p78_2_2x, eigen_quad_p78_2_3x, eigen_quad_p78_2_4x, eigen_quad_p78_2_5x, eigen_quad_p78_2_6x, eigen_quad_p78_2_7x, eigen_quad_p78_2_8x, eigen_quad_p78_3_2x, eigen_quad_p78_3_3x, eigen_quad_p78_3_4x, eigen_quad_p78_3_5x, eigen_quad_p78_3_6x, eigen_quad_p78_3_7x, eigen_quad_p78_3_8x, eigen_quad_p78_4_2x, eigen_quad_p78_4_3x, eigen_quad_p78_4_4x, eigen_quad_p78_4_5x, eigen_quad_p78_4_6x, eigen_quad_p78_4_7x, eigen_quad_p78_4_8x, eigen_quad_p78_5_2x, eigen_quad_p78_5_3x, eigen_quad_p78_5_4x, eigen_quad_p78_5_5x, eigen_quad_p78_5_6x, eigen_quad_p78_5_7x, eigen_quad_p78_5_8x, eigen_quad_p78_6_2x, eigen_quad_p78_6_3x, eigen_quad_p78_6_4x, eigen_quad_p78_6_5x, eigen_quad_p78_6_6x, eigen_quad_p78_6_7x, eigen_quad_p78_6_8x, eigen_quad_p78_7_2x, eigen_quad_p78_7_3x, eigen_quad_p78_7_4x, eigen_quad_p78_7_5x, eigen_quad_p78_7_6x, eigen_quad_p78_7_7x, eigen_quad_p78_7_8x, eigen_quad_p79_2_2x, eigen_quad_p79_2_3x, eigen_quad_p79_2_4x, eigen_quad_p79_2_5x, eigen_quad_p79_2_6x, eigen_quad_p79_2_7x, eigen_quad_p79_2_8x, eigen_quad_p79_2_9x, eigen_quad_p79_3_2x, eigen_quad_p79_3_3x, eigen_quad_p79_3_4x, eigen_quad_p79_3_5x, eigen_quad_p79_3_6x, eigen_quad_p79_3_7x, eigen_quad_p79_3_8x, eigen_quad_p79_3_9x, eigen_quad_p79_4_2x, eigen_quad_p79_4_3x, eigen_quad_p79_4_4x, eigen_quad_p79_4_5x, eigen_quad_p79_4_6x, eigen_quad_p79_4_7x, eigen_quad_p79_4_8x, eigen_quad_p79_4_9x, eigen_quad_p79_5_2x, eigen_quad_p79_5_3x, eigen_quad_p79_5_4x, eigen_quad_p79_5_5x, eigen_quad_p79_5_6x, eigen_quad_p79_5_7x, eigen_quad_p79_5_8x, eigen_quad_p79_5_9x, eigen_quad_p79_6_2x, eigen_quad_p79_6_3x, eigen_quad_p79_6_4x, eigen_quad_p79_6_5x, eigen_quad_p79_6_6x, eigen_quad_p79_6_7x, eigen_quad_p79_6_8x, eigen_quad_p79_6_9x, eigen_quad_p79_7_2x, eigen_quad_p79_7_3x, eigen_quad_p79_7_4x, eigen_quad_p79_7_5x, eigen_quad_p79_7_6x, eigen_quad_p79_7_7x, eigen_quad_p79_7_8x, eigen_quad_p79_7_9x, eigen_quad_p710_2_2x, eigen_quad_p710_2_3x, eigen_quad_p710_2_4x, eigen_quad_p710_2_5x, eigen_quad_p710_2_6x, eigen_quad_p710_2_7x, eigen_quad_p710_2_8x, eigen_quad_p710_2_9x, eigen_quad_p710_2_10x, eigen_quad_p710_3_2x, eigen_quad_p710_3_3x, eigen_quad_p710_3_4x, eigen_quad_p710_3_5x, eigen_quad_p710_3_6x, eigen_quad_p710_3_7x, eigen_quad_p710_3_8x, eigen_quad_p710_3_9x, eigen_quad_p710_3_10x, eigen_quad_p710_4_2x, eigen_quad_p710_4_3x, eigen_quad_p710_4_4x, eigen_quad_p710_4_5x, eigen_quad_p710_4_6x, eigen_quad_p710_4_7x, eigen_quad_p710_4_8x, eigen_quad_p710_4_9x, eigen_quad_p710_4_10x, eigen_quad_p710_5_2x, eigen_quad_p710_5_3x, eigen_quad_p710_5_4x, eigen_quad_p710_5_5x, eigen_quad_p710_5_6x, eigen_quad_p710_5_7x, eigen_quad_p710_5_8x, eigen_quad_p710_5_9x, eigen_quad_p710_5_10x, eigen_quad_p710_6_2x, eigen_quad_p710_6_3x, eigen_quad_p710_6_4x, eigen_quad_p710_6_5x, eigen_quad_p710_6_6x, eigen_quad_p710_6_7x, eigen_quad_p710_6_8x, eigen_quad_p710_6_9x, eigen_quad_p710_6_10x, eigen_quad_p710_7_2x, eigen_quad_p710_7_3x, eigen_quad_p710_7_4x, eigen_quad_p710_7_5x, eigen_quad_p710_7_6x, eigen_quad_p710_7_7x, eigen_quad_p710_7_8x, eigen_quad_p710_7_9x, eigen_quad_p710_7_10x,

  eigen_quad_p82_2_2x, eigen_quad_p82_3_2x, eigen_quad_p82_4_2x, eigen_quad_p82_5_2x, eigen_quad_p82_6_2x, eigen_quad_p82_7_2x, eigen_quad_p82_8_2x, eigen_quad_p83_2_2x, eigen_quad_p83_2_3x, eigen_quad_p83_3_2x, eigen_quad_p83_3_3x, eigen_quad_p83_4_2x, eigen_quad_p83_4_3x, eigen_quad_p83_5_2x, eigen_quad_p83_5_3x, eigen_quad_p83_6_2x, eigen_quad_p83_6_3x, eigen_quad_p83_7_2x, eigen_quad_p83_7_3x, eigen_quad_p83_8_2x, eigen_quad_p83_8_3x, eigen_quad_p84_2_2x, eigen_quad_p84_2_3x, eigen_quad_p84_2_4x, eigen_quad_p84_3_2x, eigen_quad_p84_3_3x, eigen_quad_p84_3_4x, eigen_quad_p84_4_2x, eigen_quad_p84_4_3x, eigen_quad_p84_4_4x, eigen_quad_p84_5_2x, eigen_quad_p84_5_3x, eigen_quad_p84_5_4x, eigen_quad_p84_6_2x, eigen_quad_p84_6_3x, eigen_quad_p84_6_4x, eigen_quad_p84_7_2x, eigen_quad_p84_7_3x, eigen_quad_p84_7_4x, eigen_quad_p84_8_2x, eigen_quad_p84_8_3x, eigen_quad_p84_8_4x, eigen_quad_p85_2_2x, eigen_quad_p85_2_3x, eigen_quad_p85_2_4x, eigen_quad_p85_2_5x, eigen_quad_p85_3_2x, eigen_quad_p85_3_3x, eigen_quad_p85_3_4x, eigen_quad_p85_3_5x, eigen_quad_p85_4_2x, eigen_quad_p85_4_3x, eigen_quad_p85_4_4x, eigen_quad_p85_4_5x, eigen_quad_p85_5_2x, eigen_quad_p85_5_3x, eigen_quad_p85_5_4x, eigen_quad_p85_5_5x, eigen_quad_p85_6_2x, eigen_quad_p85_6_3x, eigen_quad_p85_6_4x, eigen_quad_p85_6_5x, eigen_quad_p85_7_2x, eigen_quad_p85_7_3x, eigen_quad_p85_7_4x, eigen_quad_p85_7_5x, eigen_quad_p85_8_2x, eigen_quad_p85_8_3x, eigen_quad_p85_8_4x, eigen_quad_p85_8_5x, eigen_quad_p86_2_2x, eigen_quad_p86_2_3x, eigen_quad_p86_2_4x, eigen_quad_p86_2_5x, eigen_quad_p86_2_6x, eigen_quad_p86_3_2x, eigen_quad_p86_3_3x, eigen_quad_p86_3_4x, eigen_quad_p86_3_5x, eigen_quad_p86_3_6x, eigen_quad_p86_4_2x, eigen_quad_p86_4_3x, eigen_quad_p86_4_4x, eigen_quad_p86_4_5x, eigen_quad_p86_4_6x, eigen_quad_p86_5_2x, eigen_quad_p86_5_3x, eigen_quad_p86_5_4x, eigen_quad_p86_5_5x, eigen_quad_p86_5_6x, eigen_quad_p86_6_2x, eigen_quad_p86_6_3x, eigen_quad_p86_6_4x, eigen_quad_p86_6_5x, eigen_quad_p86_6_6x, eigen_quad_p86_7_2x, eigen_quad_p86_7_3x, eigen_quad_p86_7_4x, eigen_quad_p86_7_5x, eigen_quad_p86_7_6x, eigen_quad_p86_8_2x, eigen_quad_p86_8_3x, eigen_quad_p86_8_4x, eigen_quad_p86_8_5x, eigen_quad_p86_8_6x, eigen_quad_p87_2_2x, eigen_quad_p87_2_3x, eigen_quad_p87_2_4x, eigen_quad_p87_2_5x, eigen_quad_p87_2_6x, eigen_quad_p87_2_7x, eigen_quad_p87_3_2x, eigen_quad_p87_3_3x, eigen_quad_p87_3_4x, eigen_quad_p87_3_5x, eigen_quad_p87_3_6x, eigen_quad_p87_3_7x, eigen_quad_p87_4_2x, eigen_quad_p87_4_3x, eigen_quad_p87_4_4x, eigen_quad_p87_4_5x, eigen_quad_p87_4_6x, eigen_quad_p87_4_7x, eigen_quad_p87_5_2x, eigen_quad_p87_5_3x, eigen_quad_p87_5_4x, eigen_quad_p87_5_5x, eigen_quad_p87_5_6x, eigen_quad_p87_5_7x, eigen_quad_p87_6_2x, eigen_quad_p87_6_3x, eigen_quad_p87_6_4x, eigen_quad_p87_6_5x, eigen_quad_p87_6_6x, eigen_quad_p87_6_7x, eigen_quad_p87_7_2x, eigen_quad_p87_7_3x, eigen_quad_p87_7_4x, eigen_quad_p87_7_5x, eigen_quad_p87_7_6x, eigen_quad_p87_7_7x, eigen_quad_p87_8_2x, eigen_quad_p87_8_3x, eigen_quad_p87_8_4x, eigen_quad_p87_8_5x, eigen_quad_p87_8_6x, eigen_quad_p87_8_7x, eigen_quad_p88_2_2x, eigen_quad_p88_2_3x, eigen_quad_p88_2_4x, eigen_quad_p88_2_5x, eigen_quad_p88_2_6x, eigen_quad_p88_2_7x, eigen_quad_p88_2_8x, eigen_quad_p88_3_2x, eigen_quad_p88_3_3x, eigen_quad_p88_3_4x, eigen_quad_p88_3_5x, eigen_quad_p88_3_6x, eigen_quad_p88_3_7x, eigen_quad_p88_3_8x, eigen_quad_p88_4_2x, eigen_quad_p88_4_3x, eigen_quad_p88_4_4x, eigen_quad_p88_4_5x, eigen_quad_p88_4_6x, eigen_quad_p88_4_7x, eigen_quad_p88_4_8x, eigen_quad_p88_5_2x, eigen_quad_p88_5_3x, eigen_quad_p88_5_4x, eigen_quad_p88_5_5x, eigen_quad_p88_5_6x, eigen_quad_p88_5_7x, eigen_quad_p88_5_8x, eigen_quad_p88_6_2x, eigen_quad_p88_6_3x, eigen_quad_p88_6_4x, eigen_quad_p88_6_5x, eigen_quad_p88_6_6x, eigen_quad_p88_6_7x, eigen_quad_p88_6_8x, eigen_quad_p88_7_2x, eigen_quad_p88_7_3x, eigen_quad_p88_7_4x, eigen_quad_p88_7_5x, eigen_quad_p88_7_6x, eigen_quad_p88_7_7x, eigen_quad_p88_7_8x, eigen_quad_p88_8_2x, eigen_quad_p88_8_3x, eigen_quad_p88_8_4x, eigen_quad_p88_8_5x, eigen_quad_p88_8_6x, eigen_quad_p88_8_7x, eigen_quad_p88_8_8x, eigen_quad_p89_2_2x, eigen_quad_p89_2_3x, eigen_quad_p89_2_4x, eigen_quad_p89_2_5x, eigen_quad_p89_2_6x, eigen_quad_p89_2_7x, eigen_quad_p89_2_8x, eigen_quad_p89_2_9x, eigen_quad_p89_3_2x, eigen_quad_p89_3_3x, eigen_quad_p89_3_4x, eigen_quad_p89_3_5x, eigen_quad_p89_3_6x, eigen_quad_p89_3_7x, eigen_quad_p89_3_8x, eigen_quad_p89_3_9x, eigen_quad_p89_4_2x, eigen_quad_p89_4_3x, eigen_quad_p89_4_4x, eigen_quad_p89_4_5x, eigen_quad_p89_4_6x, eigen_quad_p89_4_7x, eigen_quad_p89_4_8x, eigen_quad_p89_4_9x, eigen_quad_p89_5_2x, eigen_quad_p89_5_3x, eigen_quad_p89_5_4x, eigen_quad_p89_5_5x, eigen_quad_p89_5_6x, eigen_quad_p89_5_7x, eigen_quad_p89_5_8x, eigen_quad_p89_5_9x, eigen_quad_p89_6_2x, eigen_quad_p89_6_3x, eigen_quad_p89_6_4x, eigen_quad_p89_6_5x, eigen_quad_p89_6_6x, eigen_quad_p89_6_7x, eigen_quad_p89_6_8x, eigen_quad_p89_6_9x, eigen_quad_p89_7_2x, eigen_quad_p89_7_3x, eigen_quad_p89_7_4x, eigen_quad_p89_7_5x, eigen_quad_p89_7_6x, eigen_quad_p89_7_7x, eigen_quad_p89_7_8x, eigen_quad_p89_7_9x, eigen_quad_p89_8_2x, eigen_quad_p89_8_3x, eigen_quad_p89_8_4x, eigen_quad_p89_8_5x, eigen_quad_p89_8_6x, eigen_quad_p89_8_7x, eigen_quad_p89_8_8x, eigen_quad_p89_8_9x, eigen_quad_p810_2_2x, eigen_quad_p810_2_3x, eigen_quad_p810_2_4x, eigen_quad_p810_2_5x, eigen_quad_p810_2_6x, eigen_quad_p810_2_7x, eigen_quad_p810_2_8x, eigen_quad_p810_2_9x, eigen_quad_p810_2_10x, eigen_quad_p810_3_2x, eigen_quad_p810_3_3x, eigen_quad_p810_3_4x, eigen_quad_p810_3_5x, eigen_quad_p810_3_6x, eigen_quad_p810_3_7x, eigen_quad_p810_3_8x, eigen_quad_p810_3_9x, eigen_quad_p810_3_10x, eigen_quad_p810_4_2x, eigen_quad_p810_4_3x, eigen_quad_p810_4_4x, eigen_quad_p810_4_5x, eigen_quad_p810_4_6x, eigen_quad_p810_4_7x, eigen_quad_p810_4_8x, eigen_quad_p810_4_9x, eigen_quad_p810_4_10x, eigen_quad_p810_5_2x, eigen_quad_p810_5_3x, eigen_quad_p810_5_4x, eigen_quad_p810_5_5x, eigen_quad_p810_5_6x, eigen_quad_p810_5_7x, eigen_quad_p810_5_8x, eigen_quad_p810_5_9x, eigen_quad_p810_5_10x, eigen_quad_p810_6_2x, eigen_quad_p810_6_3x, eigen_quad_p810_6_4x, eigen_quad_p810_6_5x, eigen_quad_p810_6_6x, eigen_quad_p810_6_7x, eigen_quad_p810_6_8x, eigen_quad_p810_6_9x, eigen_quad_p810_6_10x, eigen_quad_p810_7_2x, eigen_quad_p810_7_3x, eigen_quad_p810_7_4x, eigen_quad_p810_7_5x, eigen_quad_p810_7_6x, eigen_quad_p810_7_7x, eigen_quad_p810_7_8x, eigen_quad_p810_7_9x, eigen_quad_p810_7_10x, eigen_quad_p810_8_2x, eigen_quad_p810_8_3x, eigen_quad_p810_8_4x, eigen_quad_p810_8_5x, eigen_quad_p810_8_6x, eigen_quad_p810_8_7x, eigen_quad_p810_8_8x, eigen_quad_p810_8_9x, eigen_quad_p810_8_10x,

  eigen_quad_p92_2_2x, eigen_quad_p92_3_2x, eigen_quad_p92_4_2x, eigen_quad_p92_5_2x, eigen_quad_p92_6_2x, eigen_quad_p92_7_2x, eigen_quad_p92_8_2x, eigen_quad_p92_9_2x, eigen_quad_p93_2_2x, eigen_quad_p93_2_3x, eigen_quad_p93_3_2x, eigen_quad_p93_3_3x, eigen_quad_p93_4_2x, eigen_quad_p93_4_3x, eigen_quad_p93_5_2x, eigen_quad_p93_5_3x, eigen_quad_p93_6_2x, eigen_quad_p93_6_3x, eigen_quad_p93_7_2x, eigen_quad_p93_7_3x, eigen_quad_p93_8_2x, eigen_quad_p93_8_3x, eigen_quad_p93_9_2x, eigen_quad_p93_9_3x, eigen_quad_p94_2_2x, eigen_quad_p94_2_3x, eigen_quad_p94_2_4x, eigen_quad_p94_3_2x, eigen_quad_p94_3_3x, eigen_quad_p94_3_4x, eigen_quad_p94_4_2x, eigen_quad_p94_4_3x, eigen_quad_p94_4_4x, eigen_quad_p94_5_2x, eigen_quad_p94_5_3x, eigen_quad_p94_5_4x, eigen_quad_p94_6_2x, eigen_quad_p94_6_3x, eigen_quad_p94_6_4x, eigen_quad_p94_7_2x, eigen_quad_p94_7_3x, eigen_quad_p94_7_4x, eigen_quad_p94_8_2x, eigen_quad_p94_8_3x, eigen_quad_p94_8_4x, eigen_quad_p94_9_2x, eigen_quad_p94_9_3x, eigen_quad_p94_9_4x, eigen_quad_p95_2_2x, eigen_quad_p95_2_3x, eigen_quad_p95_2_4x, eigen_quad_p95_2_5x, eigen_quad_p95_3_2x, eigen_quad_p95_3_3x, eigen_quad_p95_3_4x, eigen_quad_p95_3_5x, eigen_quad_p95_4_2x, eigen_quad_p95_4_3x, eigen_quad_p95_4_4x, eigen_quad_p95_4_5x, eigen_quad_p95_5_2x, eigen_quad_p95_5_3x, eigen_quad_p95_5_4x, eigen_quad_p95_5_5x, eigen_quad_p95_6_2x, eigen_quad_p95_6_3x, eigen_quad_p95_6_4x, eigen_quad_p95_6_5x, eigen_quad_p95_7_2x, eigen_quad_p95_7_3x, eigen_quad_p95_7_4x, eigen_quad_p95_7_5x, eigen_quad_p95_8_2x, eigen_quad_p95_8_3x, eigen_quad_p95_8_4x, eigen_quad_p95_8_5x, eigen_quad_p95_9_2x, eigen_quad_p95_9_3x, eigen_quad_p95_9_4x, eigen_quad_p95_9_5x, eigen_quad_p96_2_2x, eigen_quad_p96_2_3x, eigen_quad_p96_2_4x, eigen_quad_p96_2_5x, eigen_quad_p96_2_6x, eigen_quad_p96_3_2x, eigen_quad_p96_3_3x, eigen_quad_p96_3_4x, eigen_quad_p96_3_5x, eigen_quad_p96_3_6x, eigen_quad_p96_4_2x, eigen_quad_p96_4_3x, eigen_quad_p96_4_4x, eigen_quad_p96_4_5x, eigen_quad_p96_4_6x, eigen_quad_p96_5_2x, eigen_quad_p96_5_3x, eigen_quad_p96_5_4x, eigen_quad_p96_5_5x, eigen_quad_p96_5_6x, eigen_quad_p96_6_2x, eigen_quad_p96_6_3x, eigen_quad_p96_6_4x, eigen_quad_p96_6_5x, eigen_quad_p96_6_6x, eigen_quad_p96_7_2x, eigen_quad_p96_7_3x, eigen_quad_p96_7_4x, eigen_quad_p96_7_5x, eigen_quad_p96_7_6x, eigen_quad_p96_8_2x, eigen_quad_p96_8_3x, eigen_quad_p96_8_4x, eigen_quad_p96_8_5x, eigen_quad_p96_8_6x, eigen_quad_p96_9_2x, eigen_quad_p96_9_3x, eigen_quad_p96_9_4x, eigen_quad_p96_9_5x, eigen_quad_p96_9_6x, eigen_quad_p97_2_2x, eigen_quad_p97_2_3x, eigen_quad_p97_2_4x, eigen_quad_p97_2_5x, eigen_quad_p97_2_6x, eigen_quad_p97_2_7x, eigen_quad_p97_3_2x, eigen_quad_p97_3_3x, eigen_quad_p97_3_4x, eigen_quad_p97_3_5x, eigen_quad_p97_3_6x, eigen_quad_p97_3_7x, eigen_quad_p97_4_2x, eigen_quad_p97_4_3x, eigen_quad_p97_4_4x, eigen_quad_p97_4_5x, eigen_quad_p97_4_6x, eigen_quad_p97_4_7x, eigen_quad_p97_5_2x, eigen_quad_p97_5_3x, eigen_quad_p97_5_4x, eigen_quad_p97_5_5x, eigen_quad_p97_5_6x, eigen_quad_p97_5_7x, eigen_quad_p97_6_2x, eigen_quad_p97_6_3x, eigen_quad_p97_6_4x, eigen_quad_p97_6_5x, eigen_quad_p97_6_6x, eigen_quad_p97_6_7x, eigen_quad_p97_7_2x, eigen_quad_p97_7_3x, eigen_quad_p97_7_4x, eigen_quad_p97_7_5x, eigen_quad_p97_7_6x, eigen_quad_p97_7_7x, eigen_quad_p97_8_2x, eigen_quad_p97_8_3x, eigen_quad_p97_8_4x, eigen_quad_p97_8_5x, eigen_quad_p97_8_6x, eigen_quad_p97_8_7x, eigen_quad_p97_9_2x, eigen_quad_p97_9_3x, eigen_quad_p97_9_4x, eigen_quad_p97_9_5x, eigen_quad_p97_9_6x, eigen_quad_p97_9_7x, eigen_quad_p98_2_2x, eigen_quad_p98_2_3x, eigen_quad_p98_2_4x, eigen_quad_p98_2_5x, eigen_quad_p98_2_6x, eigen_quad_p98_2_7x, eigen_quad_p98_2_8x, eigen_quad_p98_3_2x, eigen_quad_p98_3_3x, eigen_quad_p98_3_4x, eigen_quad_p98_3_5x, eigen_quad_p98_3_6x, eigen_quad_p98_3_7x, eigen_quad_p98_3_8x, eigen_quad_p98_4_2x, eigen_quad_p98_4_3x, eigen_quad_p98_4_4x, eigen_quad_p98_4_5x, eigen_quad_p98_4_6x, eigen_quad_p98_4_7x, eigen_quad_p98_4_8x, eigen_quad_p98_5_2x, eigen_quad_p98_5_3x, eigen_quad_p98_5_4x, eigen_quad_p98_5_5x, eigen_quad_p98_5_6x, eigen_quad_p98_5_7x, eigen_quad_p98_5_8x, eigen_quad_p98_6_2x, eigen_quad_p98_6_3x, eigen_quad_p98_6_4x, eigen_quad_p98_6_5x, eigen_quad_p98_6_6x, eigen_quad_p98_6_7x, eigen_quad_p98_6_8x, eigen_quad_p98_7_2x, eigen_quad_p98_7_3x, eigen_quad_p98_7_4x, eigen_quad_p98_7_5x, eigen_quad_p98_7_6x, eigen_quad_p98_7_7x, eigen_quad_p98_7_8x, eigen_quad_p98_8_2x, eigen_quad_p98_8_3x, eigen_quad_p98_8_4x, eigen_quad_p98_8_5x, eigen_quad_p98_8_6x, eigen_quad_p98_8_7x, eigen_quad_p98_8_8x, eigen_quad_p98_9_2x, eigen_quad_p98_9_3x, eigen_quad_p98_9_4x, eigen_quad_p98_9_5x, eigen_quad_p98_9_6x, eigen_quad_p98_9_7x, eigen_quad_p98_9_8x, eigen_quad_p99_2_2x, eigen_quad_p99_2_3x, eigen_quad_p99_2_4x, eigen_quad_p99_2_5x, eigen_quad_p99_2_6x, eigen_quad_p99_2_7x, eigen_quad_p99_2_8x, eigen_quad_p99_2_9x, eigen_quad_p99_3_2x, eigen_quad_p99_3_3x, eigen_quad_p99_3_4x, eigen_quad_p99_3_5x, eigen_quad_p99_3_6x, eigen_quad_p99_3_7x, eigen_quad_p99_3_8x, eigen_quad_p99_3_9x, eigen_quad_p99_4_2x, eigen_quad_p99_4_3x, eigen_quad_p99_4_4x, eigen_quad_p99_4_5x, eigen_quad_p99_4_6x, eigen_quad_p99_4_7x, eigen_quad_p99_4_8x, eigen_quad_p99_4_9x, eigen_quad_p99_5_2x, eigen_quad_p99_5_3x, eigen_quad_p99_5_4x, eigen_quad_p99_5_5x, eigen_quad_p99_5_6x, eigen_quad_p99_5_7x, eigen_quad_p99_5_8x, eigen_quad_p99_5_9x, eigen_quad_p99_6_2x, eigen_quad_p99_6_3x, eigen_quad_p99_6_4x, eigen_quad_p99_6_5x, eigen_quad_p99_6_6x, eigen_quad_p99_6_7x, eigen_quad_p99_6_8x, eigen_quad_p99_6_9x, eigen_quad_p99_7_2x, eigen_quad_p99_7_3x, eigen_quad_p99_7_4x, eigen_quad_p99_7_5x, eigen_quad_p99_7_6x, eigen_quad_p99_7_7x, eigen_quad_p99_7_8x, eigen_quad_p99_7_9x, eigen_quad_p99_8_2x, eigen_quad_p99_8_3x, eigen_quad_p99_8_4x, eigen_quad_p99_8_5x, eigen_quad_p99_8_6x, eigen_quad_p99_8_7x, eigen_quad_p99_8_8x, eigen_quad_p99_8_9x, eigen_quad_p99_9_2x, eigen_quad_p99_9_3x, eigen_quad_p99_9_4x, eigen_quad_p99_9_5x, eigen_quad_p99_9_6x, eigen_quad_p99_9_7x, eigen_quad_p99_9_8x, eigen_quad_p99_9_9x, eigen_quad_p910_2_2x, eigen_quad_p910_2_3x, eigen_quad_p910_2_4x, eigen_quad_p910_2_5x, eigen_quad_p910_2_6x, eigen_quad_p910_2_7x, eigen_quad_p910_2_8x, eigen_quad_p910_2_9x, eigen_quad_p910_2_10x, eigen_quad_p910_3_2x, eigen_quad_p910_3_3x, eigen_quad_p910_3_4x, eigen_quad_p910_3_5x, eigen_quad_p910_3_6x, eigen_quad_p910_3_7x, eigen_quad_p910_3_8x, eigen_quad_p910_3_9x, eigen_quad_p910_3_10x, eigen_quad_p910_4_2x, eigen_quad_p910_4_3x, eigen_quad_p910_4_4x, eigen_quad_p910_4_5x, eigen_quad_p910_4_6x, eigen_quad_p910_4_7x, eigen_quad_p910_4_8x, eigen_quad_p910_4_9x, eigen_quad_p910_4_10x, eigen_quad_p910_5_2x, eigen_quad_p910_5_3x, eigen_quad_p910_5_4x, eigen_quad_p910_5_5x, eigen_quad_p910_5_6x, eigen_quad_p910_5_7x, eigen_quad_p910_5_8x, eigen_quad_p910_5_9x, eigen_quad_p910_5_10x, eigen_quad_p910_6_2x, eigen_quad_p910_6_3x, eigen_quad_p910_6_4x, eigen_quad_p910_6_5x, eigen_quad_p910_6_6x, eigen_quad_p910_6_7x, eigen_quad_p910_6_8x, eigen_quad_p910_6_9x, eigen_quad_p910_6_10x, eigen_quad_p910_7_2x, eigen_quad_p910_7_3x, eigen_quad_p910_7_4x, eigen_quad_p910_7_5x, eigen_quad_p910_7_6x, eigen_quad_p910_7_7x, eigen_quad_p910_7_8x, eigen_quad_p910_7_9x, eigen_quad_p910_7_10x, eigen_quad_p910_8_2x, eigen_quad_p910_8_3x, eigen_quad_p910_8_4x, eigen_quad_p910_8_5x, eigen_quad_p910_8_6x, eigen_quad_p910_8_7x, eigen_quad_p910_8_8x, eigen_quad_p910_8_9x, eigen_quad_p910_8_10x, eigen_quad_p910_9_2x, eigen_quad_p910_9_3x, eigen_quad_p910_9_4x, eigen_quad_p910_9_5x, eigen_quad_p910_9_6x, eigen_quad_p910_9_7x, eigen_quad_p910_9_8x, eigen_quad_p910_9_9x, eigen_quad_p910_9_10x,

  eigen_quad_p102_2_2x, eigen_quad_p102_3_2x, eigen_quad_p102_4_2x, eigen_quad_p102_5_2x, eigen_quad_p102_6_2x, eigen_quad_p102_7_2x, eigen_quad_p102_8_2x, eigen_quad_p102_9_2x, eigen_quad_p102_10_2x, eigen_quad_p103_2_2x, eigen_quad_p103_2_3x, eigen_quad_p103_3_2x, eigen_quad_p103_3_3x, eigen_quad_p103_4_2x, eigen_quad_p103_4_3x, eigen_quad_p103_5_2x, eigen_quad_p103_5_3x, eigen_quad_p103_6_2x, eigen_quad_p103_6_3x, eigen_quad_p103_7_2x, eigen_quad_p103_7_3x, eigen_quad_p103_8_2x, eigen_quad_p103_8_3x, eigen_quad_p103_9_2x, eigen_quad_p103_9_3x, eigen_quad_p103_10_2x, eigen_quad_p103_10_3x, eigen_quad_p104_2_2x, eigen_quad_p104_2_3x, eigen_quad_p104_2_4x, eigen_quad_p104_3_2x, eigen_quad_p104_3_3x, eigen_quad_p104_3_4x, eigen_quad_p104_4_2x, eigen_quad_p104_4_3x, eigen_quad_p104_4_4x, eigen_quad_p104_5_2x, eigen_quad_p104_5_3x, eigen_quad_p104_5_4x, eigen_quad_p104_6_2x, eigen_quad_p104_6_3x, eigen_quad_p104_6_4x, eigen_quad_p104_7_2x, eigen_quad_p104_7_3x, eigen_quad_p104_7_4x, eigen_quad_p104_8_2x, eigen_quad_p104_8_3x, eigen_quad_p104_8_4x, eigen_quad_p104_9_2x, eigen_quad_p104_9_3x, eigen_quad_p104_9_4x, eigen_quad_p104_10_2x, eigen_quad_p104_10_3x, eigen_quad_p104_10_4x, eigen_quad_p105_2_2x, eigen_quad_p105_2_3x, eigen_quad_p105_2_4x, eigen_quad_p105_2_5x, eigen_quad_p105_3_2x, eigen_quad_p105_3_3x, eigen_quad_p105_3_4x, eigen_quad_p105_3_5x, eigen_quad_p105_4_2x, eigen_quad_p105_4_3x, eigen_quad_p105_4_4x, eigen_quad_p105_4_5x, eigen_quad_p105_5_2x, eigen_quad_p105_5_3x, eigen_quad_p105_5_4x, eigen_quad_p105_5_5x, eigen_quad_p105_6_2x, eigen_quad_p105_6_3x, eigen_quad_p105_6_4x, eigen_quad_p105_6_5x, eigen_quad_p105_7_2x, eigen_quad_p105_7_3x, eigen_quad_p105_7_4x, eigen_quad_p105_7_5x, eigen_quad_p105_8_2x, eigen_quad_p105_8_3x, eigen_quad_p105_8_4x, eigen_quad_p105_8_5x, eigen_quad_p105_9_2x, eigen_quad_p105_9_3x, eigen_quad_p105_9_4x, eigen_quad_p105_9_5x, eigen_quad_p105_10_2x, eigen_quad_p105_10_3x, eigen_quad_p105_10_4x, eigen_quad_p105_10_5x, eigen_quad_p106_2_2x, eigen_quad_p106_2_3x, eigen_quad_p106_2_4x, eigen_quad_p106_2_5x, eigen_quad_p106_2_6x, eigen_quad_p106_3_2x, eigen_quad_p106_3_3x, eigen_quad_p106_3_4x, eigen_quad_p106_3_5x, eigen_quad_p106_3_6x, eigen_quad_p106_4_2x, eigen_quad_p106_4_3x, eigen_quad_p106_4_4x, eigen_quad_p106_4_5x, eigen_quad_p106_4_6x, eigen_quad_p106_5_2x, eigen_quad_p106_5_3x, eigen_quad_p106_5_4x, eigen_quad_p106_5_5x, eigen_quad_p106_5_6x, eigen_quad_p106_6_2x, eigen_quad_p106_6_3x, eigen_quad_p106_6_4x, eigen_quad_p106_6_5x, eigen_quad_p106_6_6x, eigen_quad_p106_7_2x, eigen_quad_p106_7_3x, eigen_quad_p106_7_4x, eigen_quad_p106_7_5x, eigen_quad_p106_7_6x, eigen_quad_p106_8_2x, eigen_quad_p106_8_3x, eigen_quad_p106_8_4x, eigen_quad_p106_8_5x, eigen_quad_p106_8_6x, eigen_quad_p106_9_2x, eigen_quad_p106_9_3x, eigen_quad_p106_9_4x, eigen_quad_p106_9_5x, eigen_quad_p106_9_6x, eigen_quad_p106_10_2x, eigen_quad_p106_10_3x, eigen_quad_p106_10_4x, eigen_quad_p106_10_5x, eigen_quad_p106_10_6x, eigen_quad_p107_2_2x, eigen_quad_p107_2_3x, eigen_quad_p107_2_4x, eigen_quad_p107_2_5x, eigen_quad_p107_2_6x, eigen_quad_p107_2_7x, eigen_quad_p107_3_2x, eigen_quad_p107_3_3x, eigen_quad_p107_3_4x, eigen_quad_p107_3_5x, eigen_quad_p107_3_6x, eigen_quad_p107_3_7x, eigen_quad_p107_4_2x, eigen_quad_p107_4_3x, eigen_quad_p107_4_4x, eigen_quad_p107_4_5x, eigen_quad_p107_4_6x, eigen_quad_p107_4_7x, eigen_quad_p107_5_2x, eigen_quad_p107_5_3x, eigen_quad_p107_5_4x, eigen_quad_p107_5_5x, eigen_quad_p107_5_6x, eigen_quad_p107_5_7x, eigen_quad_p107_6_2x, eigen_quad_p107_6_3x, eigen_quad_p107_6_4x, eigen_quad_p107_6_5x, eigen_quad_p107_6_6x, eigen_quad_p107_6_7x, eigen_quad_p107_7_2x, eigen_quad_p107_7_3x, eigen_quad_p107_7_4x, eigen_quad_p107_7_5x, eigen_quad_p107_7_6x, eigen_quad_p107_7_7x, eigen_quad_p107_8_2x, eigen_quad_p107_8_3x, eigen_quad_p107_8_4x, eigen_quad_p107_8_5x, eigen_quad_p107_8_6x, eigen_quad_p107_8_7x, eigen_quad_p107_9_2x, eigen_quad_p107_9_3x, eigen_quad_p107_9_4x, eigen_quad_p107_9_5x, eigen_quad_p107_9_6x, eigen_quad_p107_9_7x, eigen_quad_p107_10_2x, eigen_quad_p107_10_3x, eigen_quad_p107_10_4x, eigen_quad_p107_10_5x, eigen_quad_p107_10_6x, eigen_quad_p107_10_7x, eigen_quad_p108_2_2x, eigen_quad_p108_2_3x, eigen_quad_p108_2_4x, eigen_quad_p108_2_5x, eigen_quad_p108_2_6x, eigen_quad_p108_2_7x, eigen_quad_p108_2_8x, eigen_quad_p108_3_2x, eigen_quad_p108_3_3x, eigen_quad_p108_3_4x, eigen_quad_p108_3_5x, eigen_quad_p108_3_6x, eigen_quad_p108_3_7x, eigen_quad_p108_3_8x, eigen_quad_p108_4_2x, eigen_quad_p108_4_3x, eigen_quad_p108_4_4x, eigen_quad_p108_4_5x, eigen_quad_p108_4_6x, eigen_quad_p108_4_7x, eigen_quad_p108_4_8x, eigen_quad_p108_5_2x, eigen_quad_p108_5_3x, eigen_quad_p108_5_4x, eigen_quad_p108_5_5x, eigen_quad_p108_5_6x, eigen_quad_p108_5_7x, eigen_quad_p108_5_8x, eigen_quad_p108_6_2x, eigen_quad_p108_6_3x, eigen_quad_p108_6_4x, eigen_quad_p108_6_5x, eigen_quad_p108_6_6x, eigen_quad_p108_6_7x, eigen_quad_p108_6_8x, eigen_quad_p108_7_2x, eigen_quad_p108_7_3x, eigen_quad_p108_7_4x, eigen_quad_p108_7_5x, eigen_quad_p108_7_6x, eigen_quad_p108_7_7x, eigen_quad_p108_7_8x, eigen_quad_p108_8_2x, eigen_quad_p108_8_3x, eigen_quad_p108_8_4x, eigen_quad_p108_8_5x, eigen_quad_p108_8_6x, eigen_quad_p108_8_7x, eigen_quad_p108_8_8x, eigen_quad_p108_9_2x, eigen_quad_p108_9_3x, eigen_quad_p108_9_4x, eigen_quad_p108_9_5x, eigen_quad_p108_9_6x, eigen_quad_p108_9_7x, eigen_quad_p108_9_8x, eigen_quad_p108_10_2x, eigen_quad_p108_10_3x, eigen_quad_p108_10_4x, eigen_quad_p108_10_5x, eigen_quad_p108_10_6x, eigen_quad_p108_10_7x, eigen_quad_p108_10_8x, eigen_quad_p109_2_2x, eigen_quad_p109_2_3x, eigen_quad_p109_2_4x, eigen_quad_p109_2_5x, eigen_quad_p109_2_6x, eigen_quad_p109_2_7x, eigen_quad_p109_2_8x, eigen_quad_p109_2_9x, eigen_quad_p109_3_2x, eigen_quad_p109_3_3x, eigen_quad_p109_3_4x, eigen_quad_p109_3_5x, eigen_quad_p109_3_6x, eigen_quad_p109_3_7x, eigen_quad_p109_3_8x, eigen_quad_p109_3_9x, eigen_quad_p109_4_2x, eigen_quad_p109_4_3x, eigen_quad_p109_4_4x, eigen_quad_p109_4_5x, eigen_quad_p109_4_6x, eigen_quad_p109_4_7x, eigen_quad_p109_4_8x, eigen_quad_p109_4_9x, eigen_quad_p109_5_2x, eigen_quad_p109_5_3x, eigen_quad_p109_5_4x, eigen_quad_p109_5_5x, eigen_quad_p109_5_6x, eigen_quad_p109_5_7x, eigen_quad_p109_5_8x, eigen_quad_p109_5_9x, eigen_quad_p109_6_2x, eigen_quad_p109_6_3x, eigen_quad_p109_6_4x, eigen_quad_p109_6_5x, eigen_quad_p109_6_6x, eigen_quad_p109_6_7x, eigen_quad_p109_6_8x, eigen_quad_p109_6_9x, eigen_quad_p109_7_2x, eigen_quad_p109_7_3x, eigen_quad_p109_7_4x, eigen_quad_p109_7_5x, eigen_quad_p109_7_6x, eigen_quad_p109_7_7x, eigen_quad_p109_7_8x, eigen_quad_p109_7_9x, eigen_quad_p109_8_2x, eigen_quad_p109_8_3x, eigen_quad_p109_8_4x, eigen_quad_p109_8_5x, eigen_quad_p109_8_6x, eigen_quad_p109_8_7x, eigen_quad_p109_8_8x, eigen_quad_p109_8_9x, eigen_quad_p109_9_2x, eigen_quad_p109_9_3x, eigen_quad_p109_9_4x, eigen_quad_p109_9_5x, eigen_quad_p109_9_6x, eigen_quad_p109_9_7x, eigen_quad_p109_9_8x, eigen_quad_p109_9_9x, eigen_quad_p109_10_2x, eigen_quad_p109_10_3x, eigen_quad_p109_10_4x, eigen_quad_p109_10_5x, eigen_quad_p109_10_6x, eigen_quad_p109_10_7x, eigen_quad_p109_10_8x, eigen_quad_p109_10_9x, eigen_quad_p1010_2_2x, eigen_quad_p1010_2_3x, eigen_quad_p1010_2_4x, eigen_quad_p1010_2_5x, eigen_quad_p1010_2_6x, eigen_quad_p1010_2_7x, eigen_quad_p1010_2_8x, eigen_quad_p1010_2_9x, eigen_quad_p1010_2_10x, eigen_quad_p1010_3_2x, eigen_quad_p1010_3_3x, eigen_quad_p1010_3_4x, eigen_quad_p1010_3_5x, eigen_quad_p1010_3_6x, eigen_quad_p1010_3_7x, eigen_quad_p1010_3_8x, eigen_quad_p1010_3_9x, eigen_quad_p1010_3_10x, eigen_quad_p1010_4_2x, eigen_quad_p1010_4_3x, eigen_quad_p1010_4_4x, eigen_quad_p1010_4_5x, eigen_quad_p1010_4_6x, eigen_quad_p1010_4_7x, eigen_quad_p1010_4_8x, eigen_quad_p1010_4_9x, eigen_quad_p1010_4_10x, eigen_quad_p1010_5_2x, eigen_quad_p1010_5_3x, eigen_quad_p1010_5_4x, eigen_quad_p1010_5_5x, eigen_quad_p1010_5_6x, eigen_quad_p1010_5_7x, eigen_quad_p1010_5_8x, eigen_quad_p1010_5_9x, eigen_quad_p1010_5_10x, eigen_quad_p1010_6_2x, eigen_quad_p1010_6_3x, eigen_quad_p1010_6_4x, eigen_quad_p1010_6_5x, eigen_quad_p1010_6_6x, eigen_quad_p1010_6_7x, eigen_quad_p1010_6_8x, eigen_quad_p1010_6_9x, eigen_quad_p1010_6_10x, eigen_quad_p1010_7_2x, eigen_quad_p1010_7_3x, eigen_quad_p1010_7_4x, eigen_quad_p1010_7_5x, eigen_quad_p1010_7_6x, eigen_quad_p1010_7_7x, eigen_quad_p1010_7_8x, eigen_quad_p1010_7_9x, eigen_quad_p1010_7_10x, eigen_quad_p1010_8_2x, eigen_quad_p1010_8_3x, eigen_quad_p1010_8_4x, eigen_quad_p1010_8_5x, eigen_quad_p1010_8_6x, eigen_quad_p1010_8_7x, eigen_quad_p1010_8_8x, eigen_quad_p1010_8_9x, eigen_quad_p1010_8_10x, eigen_quad_p1010_9_2x, eigen_quad_p1010_9_3x, eigen_quad_p1010_9_4x, eigen_quad_p1010_9_5x, eigen_quad_p1010_9_6x, eigen_quad_p1010_9_7x, eigen_quad_p1010_9_8x, eigen_quad_p1010_9_9x, eigen_quad_p1010_9_10x, eigen_quad_p1010_10_2x, eigen_quad_p1010_10_3x, eigen_quad_p1010_10_4x, eigen_quad_p1010_10_5x, eigen_quad_p1010_10_6x, eigen_quad_p1010_10_7x, eigen_quad_p1010_10_8x, eigen_quad_p1010_10_9x, eigen_quad_p1010_10_10x
};

static Shapeset::shape_fn_t eigen_quad_fn_dy[] =
{
eigen_quad_0_0y,   eigen_quad_0_1y,   eigen_quad_0_2y,   eigen_quad_0_3y_0, eigen_quad_0_3y_1,  eigen_quad_0_4y,   eigen_quad_0_5y_0,  eigen_quad_0_5y_1, eigen_quad_0_6y,   eigen_quad_0_7y_0, eigen_quad_0_7y_1, eigen_quad_0_8y,  eigen_quad_0_9y_0,  eigen_quad_0_9y_1, eigen_quad_0_10y,  eigen_quad_1_0y,  eigen_quad_1_1y,   eigen_quad_1_2y,   eigen_quad_1_3y_0, eigen_quad_1_3y_1, eigen_quad_1_4y,    eigen_quad_1_5y_0, eigen_quad_1_5y_1, eigen_quad_1_6y,  eigen_quad_1_7y_0, eigen_quad_1_7y_1, eigen_quad_1_8y,   eigen_quad_1_9y_0, eigen_quad_1_9y_1,  eigen_quad_1_10y,  eigen_quad_2_0y,   eigen_quad_2_1y,   eigen_quad_3_0y_0, eigen_quad_3_0y_1, eigen_quad_3_1y_0, eigen_quad_3_1y_1, eigen_quad_4_0y,   eigen_quad_4_1y,   eigen_quad_5_0y_0, eigen_quad_5_0y_1, eigen_quad_5_1y_0, eigen_quad_5_1y_1, eigen_quad_6_0y,   eigen_quad_6_1y,   eigen_quad_7_0y_0, eigen_quad_7_0y_1, eigen_quad_7_1y_0, eigen_quad_7_1y_1, eigen_quad_8_0y,   eigen_quad_8_1y,   eigen_quad_9_0y_0, eigen_quad_9_0y_1, eigen_quad_9_1y_0, eigen_quad_9_1y_1, eigen_quad_10_0y,  eigen_quad_10_1y,

  eigen_quad_p22_2_2y, eigen_quad_p23_2_2y, eigen_quad_p23_2_3y, eigen_quad_p24_2_2y, eigen_quad_p24_2_3y, eigen_quad_p24_2_4y, eigen_quad_p25_2_2y, eigen_quad_p25_2_3y, eigen_quad_p25_2_4y, eigen_quad_p25_2_5y, eigen_quad_p26_2_2y, eigen_quad_p26_2_3y, eigen_quad_p26_2_4y, eigen_quad_p26_2_5y, eigen_quad_p26_2_6y, eigen_quad_p27_2_2y, eigen_quad_p27_2_3y, eigen_quad_p27_2_4y, eigen_quad_p27_2_5y, eigen_quad_p27_2_6y, eigen_quad_p27_2_7y, eigen_quad_p28_2_2y, eigen_quad_p28_2_3y, eigen_quad_p28_2_4y, eigen_quad_p28_2_5y, eigen_quad_p28_2_6y, eigen_quad_p28_2_7y, eigen_quad_p28_2_8y, eigen_quad_p29_2_2y, eigen_quad_p29_2_3y, eigen_quad_p29_2_4y, eigen_quad_p29_2_5y, eigen_quad_p29_2_6y, eigen_quad_p29_2_7y, eigen_quad_p29_2_8y, eigen_quad_p29_2_9y, eigen_quad_p210_2_2y, eigen_quad_p210_2_3y, eigen_quad_p210_2_4y, eigen_quad_p210_2_5y, eigen_quad_p210_2_6y, eigen_quad_p210_2_7y, eigen_quad_p210_2_8y, eigen_quad_p210_2_9y, eigen_quad_p210_2_10y,

  eigen_quad_p32_2_2y, eigen_quad_p32_3_2y, eigen_quad_p33_2_2y, eigen_quad_p33_2_3y, eigen_quad_p33_3_2y, eigen_quad_p33_3_3y, eigen_quad_p34_2_2y, eigen_quad_p34_2_3y, eigen_quad_p34_2_4y, eigen_quad_p34_3_2y, eigen_quad_p34_3_3y, eigen_quad_p34_3_4y, eigen_quad_p35_2_2y, eigen_quad_p35_2_3y, eigen_quad_p35_2_4y, eigen_quad_p35_2_5y, eigen_quad_p35_3_2y, eigen_quad_p35_3_3y, eigen_quad_p35_3_4y, eigen_quad_p35_3_5y, eigen_quad_p36_2_2y, eigen_quad_p36_2_3y, eigen_quad_p36_2_4y, eigen_quad_p36_2_5y, eigen_quad_p36_2_6y, eigen_quad_p36_3_2y, eigen_quad_p36_3_3y, eigen_quad_p36_3_4y, eigen_quad_p36_3_5y, eigen_quad_p36_3_6y, eigen_quad_p37_2_2y, eigen_quad_p37_2_3y, eigen_quad_p37_2_4y, eigen_quad_p37_2_5y, eigen_quad_p37_2_6y, eigen_quad_p37_2_7y, eigen_quad_p37_3_2y, eigen_quad_p37_3_3y, eigen_quad_p37_3_4y, eigen_quad_p37_3_5y, eigen_quad_p37_3_6y, eigen_quad_p37_3_7y, eigen_quad_p38_2_2y, eigen_quad_p38_2_3y, eigen_quad_p38_2_4y, eigen_quad_p38_2_5y, eigen_quad_p38_2_6y, eigen_quad_p38_2_7y, eigen_quad_p38_2_8y, eigen_quad_p38_3_2y, eigen_quad_p38_3_3y, eigen_quad_p38_3_4y, eigen_quad_p38_3_5y, eigen_quad_p38_3_6y, eigen_quad_p38_3_7y, eigen_quad_p38_3_8y, eigen_quad_p39_2_2y, eigen_quad_p39_2_3y, eigen_quad_p39_2_4y, eigen_quad_p39_2_5y, eigen_quad_p39_2_6y, eigen_quad_p39_2_7y, eigen_quad_p39_2_8y, eigen_quad_p39_2_9y, eigen_quad_p39_3_2y, eigen_quad_p39_3_3y, eigen_quad_p39_3_4y, eigen_quad_p39_3_5y, eigen_quad_p39_3_6y, eigen_quad_p39_3_7y, eigen_quad_p39_3_8y, eigen_quad_p39_3_9y, eigen_quad_p310_2_2y, eigen_quad_p310_2_3y, eigen_quad_p310_2_4y, eigen_quad_p310_2_5y, eigen_quad_p310_2_6y, eigen_quad_p310_2_7y, eigen_quad_p310_2_8y, eigen_quad_p310_2_9y, eigen_quad_p310_2_10y, eigen_quad_p310_3_2y, eigen_quad_p310_3_3y, eigen_quad_p310_3_4y, eigen_quad_p310_3_5y, eigen_quad_p310_3_6y, eigen_quad_p310_3_7y, eigen_quad_p310_3_8y, eigen_quad_p310_3_9y, eigen_quad_p310_3_10y,

  eigen_quad_p42_2_2y, eigen_quad_p42_3_2y, eigen_quad_p42_4_2y, eigen_quad_p43_2_2y, eigen_quad_p43_2_3y, eigen_quad_p43_3_2y, eigen_quad_p43_3_3y, eigen_quad_p43_4_2y, eigen_quad_p43_4_3y, eigen_quad_p44_2_2y, eigen_quad_p44_2_3y, eigen_quad_p44_2_4y, eigen_quad_p44_3_2y, eigen_quad_p44_3_3y, eigen_quad_p44_3_4y, eigen_quad_p44_4_2y, eigen_quad_p44_4_3y, eigen_quad_p44_4_4y, eigen_quad_p45_2_2y, eigen_quad_p45_2_3y, eigen_quad_p45_2_4y, eigen_quad_p45_2_5y, eigen_quad_p45_3_2y, eigen_quad_p45_3_3y, eigen_quad_p45_3_4y, eigen_quad_p45_3_5y, eigen_quad_p45_4_2y, eigen_quad_p45_4_3y, eigen_quad_p45_4_4y, eigen_quad_p45_4_5y, eigen_quad_p46_2_2y, eigen_quad_p46_2_3y, eigen_quad_p46_2_4y, eigen_quad_p46_2_5y, eigen_quad_p46_2_6y, eigen_quad_p46_3_2y, eigen_quad_p46_3_3y, eigen_quad_p46_3_4y, eigen_quad_p46_3_5y, eigen_quad_p46_3_6y, eigen_quad_p46_4_2y, eigen_quad_p46_4_3y, eigen_quad_p46_4_4y, eigen_quad_p46_4_5y, eigen_quad_p46_4_6y, eigen_quad_p47_2_2y, eigen_quad_p47_2_3y, eigen_quad_p47_2_4y, eigen_quad_p47_2_5y, eigen_quad_p47_2_6y, eigen_quad_p47_2_7y, eigen_quad_p47_3_2y, eigen_quad_p47_3_3y, eigen_quad_p47_3_4y, eigen_quad_p47_3_5y, eigen_quad_p47_3_6y, eigen_quad_p47_3_7y, eigen_quad_p47_4_2y, eigen_quad_p47_4_3y, eigen_quad_p47_4_4y, eigen_quad_p47_4_5y, eigen_quad_p47_4_6y, eigen_quad_p47_4_7y, eigen_quad_p48_2_2y, eigen_quad_p48_2_3y, eigen_quad_p48_2_4y, eigen_quad_p48_2_5y, eigen_quad_p48_2_6y, eigen_quad_p48_2_7y, eigen_quad_p48_2_8y, eigen_quad_p48_3_2y, eigen_quad_p48_3_3y, eigen_quad_p48_3_4y, eigen_quad_p48_3_5y, eigen_quad_p48_3_6y, eigen_quad_p48_3_7y, eigen_quad_p48_3_8y, eigen_quad_p48_4_2y, eigen_quad_p48_4_3y, eigen_quad_p48_4_4y, eigen_quad_p48_4_5y, eigen_quad_p48_4_6y, eigen_quad_p48_4_7y, eigen_quad_p48_4_8y, eigen_quad_p49_2_2y, eigen_quad_p49_2_3y, eigen_quad_p49_2_4y, eigen_quad_p49_2_5y, eigen_quad_p49_2_6y, eigen_quad_p49_2_7y, eigen_quad_p49_2_8y, eigen_quad_p49_2_9y, eigen_quad_p49_3_2y, eigen_quad_p49_3_3y, eigen_quad_p49_3_4y, eigen_quad_p49_3_5y, eigen_quad_p49_3_6y, eigen_quad_p49_3_7y, eigen_quad_p49_3_8y, eigen_quad_p49_3_9y, eigen_quad_p49_4_2y, eigen_quad_p49_4_3y, eigen_quad_p49_4_4y, eigen_quad_p49_4_5y, eigen_quad_p49_4_6y, eigen_quad_p49_4_7y, eigen_quad_p49_4_8y, eigen_quad_p49_4_9y, eigen_quad_p410_2_2y, eigen_quad_p410_2_3y, eigen_quad_p410_2_4y, eigen_quad_p410_2_5y, eigen_quad_p410_2_6y, eigen_quad_p410_2_7y, eigen_quad_p410_2_8y, eigen_quad_p410_2_9y, eigen_quad_p410_2_10y, eigen_quad_p410_3_2y, eigen_quad_p410_3_3y, eigen_quad_p410_3_4y, eigen_quad_p410_3_5y, eigen_quad_p410_3_6y, eigen_quad_p410_3_7y, eigen_quad_p410_3_8y, eigen_quad_p410_3_9y, eigen_quad_p410_3_10y, eigen_quad_p410_4_2y, eigen_quad_p410_4_3y, eigen_quad_p410_4_4y, eigen_quad_p410_4_5y, eigen_quad_p410_4_6y, eigen_quad_p410_4_7y, eigen_quad_p410_4_8y, eigen_quad_p410_4_9y, eigen_quad_p410_4_10y,

  eigen_quad_p52_2_2y, eigen_quad_p52_3_2y, eigen_quad_p52_4_2y, eigen_quad_p52_5_2y, eigen_quad_p53_2_2y, eigen_quad_p53_2_3y, eigen_quad_p53_3_2y, eigen_quad_p53_3_3y, eigen_quad_p53_4_2y, eigen_quad_p53_4_3y, eigen_quad_p53_5_2y, eigen_quad_p53_5_3y, eigen_quad_p54_2_2y, eigen_quad_p54_2_3y, eigen_quad_p54_2_4y, eigen_quad_p54_3_2y, eigen_quad_p54_3_3y, eigen_quad_p54_3_4y, eigen_quad_p54_4_2y, eigen_quad_p54_4_3y, eigen_quad_p54_4_4y, eigen_quad_p54_5_2y, eigen_quad_p54_5_3y, eigen_quad_p54_5_4y, eigen_quad_p55_2_2y, eigen_quad_p55_2_3y, eigen_quad_p55_2_4y, eigen_quad_p55_2_5y, eigen_quad_p55_3_2y, eigen_quad_p55_3_3y, eigen_quad_p55_3_4y, eigen_quad_p55_3_5y, eigen_quad_p55_4_2y, eigen_quad_p55_4_3y, eigen_quad_p55_4_4y, eigen_quad_p55_4_5y, eigen_quad_p55_5_2y, eigen_quad_p55_5_3y, eigen_quad_p55_5_4y, eigen_quad_p55_5_5y, eigen_quad_p56_2_2y, eigen_quad_p56_2_3y, eigen_quad_p56_2_4y, eigen_quad_p56_2_5y, eigen_quad_p56_2_6y, eigen_quad_p56_3_2y, eigen_quad_p56_3_3y, eigen_quad_p56_3_4y, eigen_quad_p56_3_5y, eigen_quad_p56_3_6y, eigen_quad_p56_4_2y, eigen_quad_p56_4_3y, eigen_quad_p56_4_4y, eigen_quad_p56_4_5y, eigen_quad_p56_4_6y, eigen_quad_p56_5_2y, eigen_quad_p56_5_3y, eigen_quad_p56_5_4y, eigen_quad_p56_5_5y, eigen_quad_p56_5_6y, eigen_quad_p57_2_2y, eigen_quad_p57_2_3y, eigen_quad_p57_2_4y, eigen_quad_p57_2_5y, eigen_quad_p57_2_6y, eigen_quad_p57_2_7y, eigen_quad_p57_3_2y, eigen_quad_p57_3_3y, eigen_quad_p57_3_4y, eigen_quad_p57_3_5y, eigen_quad_p57_3_6y, eigen_quad_p57_3_7y, eigen_quad_p57_4_2y, eigen_quad_p57_4_3y, eigen_quad_p57_4_4y, eigen_quad_p57_4_5y, eigen_quad_p57_4_6y, eigen_quad_p57_4_7y, eigen_quad_p57_5_2y, eigen_quad_p57_5_3y, eigen_quad_p57_5_4y, eigen_quad_p57_5_5y, eigen_quad_p57_5_6y, eigen_quad_p57_5_7y, eigen_quad_p58_2_2y, eigen_quad_p58_2_3y, eigen_quad_p58_2_4y, eigen_quad_p58_2_5y, eigen_quad_p58_2_6y, eigen_quad_p58_2_7y, eigen_quad_p58_2_8y, eigen_quad_p58_3_2y, eigen_quad_p58_3_3y, eigen_quad_p58_3_4y, eigen_quad_p58_3_5y, eigen_quad_p58_3_6y, eigen_quad_p58_3_7y, eigen_quad_p58_3_8y, eigen_quad_p58_4_2y, eigen_quad_p58_4_3y, eigen_quad_p58_4_4y, eigen_quad_p58_4_5y, eigen_quad_p58_4_6y, eigen_quad_p58_4_7y, eigen_quad_p58_4_8y, eigen_quad_p58_5_2y, eigen_quad_p58_5_3y, eigen_quad_p58_5_4y, eigen_quad_p58_5_5y, eigen_quad_p58_5_6y, eigen_quad_p58_5_7y, eigen_quad_p58_5_8y, eigen_quad_p59_2_2y, eigen_quad_p59_2_3y, eigen_quad_p59_2_4y, eigen_quad_p59_2_5y, eigen_quad_p59_2_6y, eigen_quad_p59_2_7y, eigen_quad_p59_2_8y, eigen_quad_p59_2_9y, eigen_quad_p59_3_2y, eigen_quad_p59_3_3y, eigen_quad_p59_3_4y, eigen_quad_p59_3_5y, eigen_quad_p59_3_6y, eigen_quad_p59_3_7y, eigen_quad_p59_3_8y, eigen_quad_p59_3_9y, eigen_quad_p59_4_2y, eigen_quad_p59_4_3y, eigen_quad_p59_4_4y, eigen_quad_p59_4_5y, eigen_quad_p59_4_6y, eigen_quad_p59_4_7y, eigen_quad_p59_4_8y, eigen_quad_p59_4_9y, eigen_quad_p59_5_2y, eigen_quad_p59_5_3y, eigen_quad_p59_5_4y, eigen_quad_p59_5_5y, eigen_quad_p59_5_6y, eigen_quad_p59_5_7y, eigen_quad_p59_5_8y, eigen_quad_p59_5_9y, eigen_quad_p510_2_2y, eigen_quad_p510_2_3y, eigen_quad_p510_2_4y, eigen_quad_p510_2_5y, eigen_quad_p510_2_6y, eigen_quad_p510_2_7y, eigen_quad_p510_2_8y, eigen_quad_p510_2_9y, eigen_quad_p510_2_10y, eigen_quad_p510_3_2y, eigen_quad_p510_3_3y, eigen_quad_p510_3_4y, eigen_quad_p510_3_5y, eigen_quad_p510_3_6y, eigen_quad_p510_3_7y, eigen_quad_p510_3_8y, eigen_quad_p510_3_9y, eigen_quad_p510_3_10y, eigen_quad_p510_4_2y, eigen_quad_p510_4_3y, eigen_quad_p510_4_4y, eigen_quad_p510_4_5y, eigen_quad_p510_4_6y, eigen_quad_p510_4_7y, eigen_quad_p510_4_8y, eigen_quad_p510_4_9y, eigen_quad_p510_4_10y, eigen_quad_p510_5_2y, eigen_quad_p510_5_3y, eigen_quad_p510_5_4y, eigen_quad_p510_5_5y, eigen_quad_p510_5_6y, eigen_quad_p510_5_7y, eigen_quad_p510_5_8y, eigen_quad_p510_5_9y, eigen_quad_p510_5_10y,

  eigen_quad_p62_2_2y, eigen_quad_p62_3_2y, eigen_quad_p62_4_2y, eigen_quad_p62_5_2y, eigen_quad_p62_6_2y, eigen_quad_p63_2_2y, eigen_quad_p63_2_3y, eigen_quad_p63_3_2y, eigen_quad_p63_3_3y, eigen_quad_p63_4_2y, eigen_quad_p63_4_3y, eigen_quad_p63_5_2y, eigen_quad_p63_5_3y, eigen_quad_p63_6_2y, eigen_quad_p63_6_3y, eigen_quad_p64_2_2y, eigen_quad_p64_2_3y, eigen_quad_p64_2_4y, eigen_quad_p64_3_2y, eigen_quad_p64_3_3y, eigen_quad_p64_3_4y, eigen_quad_p64_4_2y, eigen_quad_p64_4_3y, eigen_quad_p64_4_4y, eigen_quad_p64_5_2y, eigen_quad_p64_5_3y, eigen_quad_p64_5_4y, eigen_quad_p64_6_2y, eigen_quad_p64_6_3y, eigen_quad_p64_6_4y, eigen_quad_p65_2_2y, eigen_quad_p65_2_3y, eigen_quad_p65_2_4y, eigen_quad_p65_2_5y, eigen_quad_p65_3_2y, eigen_quad_p65_3_3y, eigen_quad_p65_3_4y, eigen_quad_p65_3_5y, eigen_quad_p65_4_2y, eigen_quad_p65_4_3y, eigen_quad_p65_4_4y, eigen_quad_p65_4_5y, eigen_quad_p65_5_2y, eigen_quad_p65_5_3y, eigen_quad_p65_5_4y, eigen_quad_p65_5_5y, eigen_quad_p65_6_2y, eigen_quad_p65_6_3y, eigen_quad_p65_6_4y, eigen_quad_p65_6_5y, eigen_quad_p66_2_2y, eigen_quad_p66_2_3y, eigen_quad_p66_2_4y, eigen_quad_p66_2_5y, eigen_quad_p66_2_6y, eigen_quad_p66_3_2y, eigen_quad_p66_3_3y, eigen_quad_p66_3_4y, eigen_quad_p66_3_5y, eigen_quad_p66_3_6y, eigen_quad_p66_4_2y, eigen_quad_p66_4_3y, eigen_quad_p66_4_4y, eigen_quad_p66_4_5y, eigen_quad_p66_4_6y, eigen_quad_p66_5_2y, eigen_quad_p66_5_3y, eigen_quad_p66_5_4y, eigen_quad_p66_5_5y, eigen_quad_p66_5_6y, eigen_quad_p66_6_2y, eigen_quad_p66_6_3y, eigen_quad_p66_6_4y, eigen_quad_p66_6_5y, eigen_quad_p66_6_6y, eigen_quad_p67_2_2y, eigen_quad_p67_2_3y, eigen_quad_p67_2_4y, eigen_quad_p67_2_5y, eigen_quad_p67_2_6y, eigen_quad_p67_2_7y, eigen_quad_p67_3_2y, eigen_quad_p67_3_3y, eigen_quad_p67_3_4y, eigen_quad_p67_3_5y, eigen_quad_p67_3_6y, eigen_quad_p67_3_7y, eigen_quad_p67_4_2y, eigen_quad_p67_4_3y, eigen_quad_p67_4_4y, eigen_quad_p67_4_5y, eigen_quad_p67_4_6y, eigen_quad_p67_4_7y, eigen_quad_p67_5_2y, eigen_quad_p67_5_3y, eigen_quad_p67_5_4y, eigen_quad_p67_5_5y, eigen_quad_p67_5_6y, eigen_quad_p67_5_7y, eigen_quad_p67_6_2y, eigen_quad_p67_6_3y, eigen_quad_p67_6_4y, eigen_quad_p67_6_5y, eigen_quad_p67_6_6y, eigen_quad_p67_6_7y, eigen_quad_p68_2_2y, eigen_quad_p68_2_3y, eigen_quad_p68_2_4y, eigen_quad_p68_2_5y, eigen_quad_p68_2_6y, eigen_quad_p68_2_7y, eigen_quad_p68_2_8y, eigen_quad_p68_3_2y, eigen_quad_p68_3_3y, eigen_quad_p68_3_4y, eigen_quad_p68_3_5y, eigen_quad_p68_3_6y, eigen_quad_p68_3_7y, eigen_quad_p68_3_8y, eigen_quad_p68_4_2y, eigen_quad_p68_4_3y, eigen_quad_p68_4_4y, eigen_quad_p68_4_5y, eigen_quad_p68_4_6y, eigen_quad_p68_4_7y, eigen_quad_p68_4_8y, eigen_quad_p68_5_2y, eigen_quad_p68_5_3y, eigen_quad_p68_5_4y, eigen_quad_p68_5_5y, eigen_quad_p68_5_6y, eigen_quad_p68_5_7y, eigen_quad_p68_5_8y, eigen_quad_p68_6_2y, eigen_quad_p68_6_3y, eigen_quad_p68_6_4y, eigen_quad_p68_6_5y, eigen_quad_p68_6_6y, eigen_quad_p68_6_7y, eigen_quad_p68_6_8y, eigen_quad_p69_2_2y, eigen_quad_p69_2_3y, eigen_quad_p69_2_4y, eigen_quad_p69_2_5y, eigen_quad_p69_2_6y, eigen_quad_p69_2_7y, eigen_quad_p69_2_8y, eigen_quad_p69_2_9y, eigen_quad_p69_3_2y, eigen_quad_p69_3_3y, eigen_quad_p69_3_4y, eigen_quad_p69_3_5y, eigen_quad_p69_3_6y, eigen_quad_p69_3_7y, eigen_quad_p69_3_8y, eigen_quad_p69_3_9y, eigen_quad_p69_4_2y, eigen_quad_p69_4_3y, eigen_quad_p69_4_4y, eigen_quad_p69_4_5y, eigen_quad_p69_4_6y, eigen_quad_p69_4_7y, eigen_quad_p69_4_8y, eigen_quad_p69_4_9y, eigen_quad_p69_5_2y, eigen_quad_p69_5_3y, eigen_quad_p69_5_4y, eigen_quad_p69_5_5y, eigen_quad_p69_5_6y, eigen_quad_p69_5_7y, eigen_quad_p69_5_8y, eigen_quad_p69_5_9y, eigen_quad_p69_6_2y, eigen_quad_p69_6_3y, eigen_quad_p69_6_4y, eigen_quad_p69_6_5y, eigen_quad_p69_6_6y, eigen_quad_p69_6_7y, eigen_quad_p69_6_8y, eigen_quad_p69_6_9y, eigen_quad_p610_2_2y, eigen_quad_p610_2_3y, eigen_quad_p610_2_4y, eigen_quad_p610_2_5y, eigen_quad_p610_2_6y, eigen_quad_p610_2_7y, eigen_quad_p610_2_8y, eigen_quad_p610_2_9y, eigen_quad_p610_2_10y, eigen_quad_p610_3_2y, eigen_quad_p610_3_3y, eigen_quad_p610_3_4y, eigen_quad_p610_3_5y, eigen_quad_p610_3_6y, eigen_quad_p610_3_7y, eigen_quad_p610_3_8y, eigen_quad_p610_3_9y, eigen_quad_p610_3_10y, eigen_quad_p610_4_2y, eigen_quad_p610_4_3y, eigen_quad_p610_4_4y, eigen_quad_p610_4_5y, eigen_quad_p610_4_6y, eigen_quad_p610_4_7y, eigen_quad_p610_4_8y, eigen_quad_p610_4_9y, eigen_quad_p610_4_10y, eigen_quad_p610_5_2y, eigen_quad_p610_5_3y, eigen_quad_p610_5_4y, eigen_quad_p610_5_5y, eigen_quad_p610_5_6y, eigen_quad_p610_5_7y, eigen_quad_p610_5_8y, eigen_quad_p610_5_9y, eigen_quad_p610_5_10y, eigen_quad_p610_6_2y, eigen_quad_p610_6_3y, eigen_quad_p610_6_4y, eigen_quad_p610_6_5y, eigen_quad_p610_6_6y, eigen_quad_p610_6_7y, eigen_quad_p610_6_8y, eigen_quad_p610_6_9y, eigen_quad_p610_6_10y,

  eigen_quad_p72_2_2y, eigen_quad_p72_3_2y, eigen_quad_p72_4_2y, eigen_quad_p72_5_2y, eigen_quad_p72_6_2y, eigen_quad_p72_7_2y, eigen_quad_p73_2_2y, eigen_quad_p73_2_3y, eigen_quad_p73_3_2y, eigen_quad_p73_3_3y, eigen_quad_p73_4_2y, eigen_quad_p73_4_3y, eigen_quad_p73_5_2y, eigen_quad_p73_5_3y, eigen_quad_p73_6_2y, eigen_quad_p73_6_3y, eigen_quad_p73_7_2y, eigen_quad_p73_7_3y, eigen_quad_p74_2_2y, eigen_quad_p74_2_3y, eigen_quad_p74_2_4y, eigen_quad_p74_3_2y, eigen_quad_p74_3_3y, eigen_quad_p74_3_4y, eigen_quad_p74_4_2y, eigen_quad_p74_4_3y, eigen_quad_p74_4_4y, eigen_quad_p74_5_2y, eigen_quad_p74_5_3y, eigen_quad_p74_5_4y, eigen_quad_p74_6_2y, eigen_quad_p74_6_3y, eigen_quad_p74_6_4y, eigen_quad_p74_7_2y, eigen_quad_p74_7_3y, eigen_quad_p74_7_4y, eigen_quad_p75_2_2y, eigen_quad_p75_2_3y, eigen_quad_p75_2_4y, eigen_quad_p75_2_5y, eigen_quad_p75_3_2y, eigen_quad_p75_3_3y, eigen_quad_p75_3_4y, eigen_quad_p75_3_5y, eigen_quad_p75_4_2y, eigen_quad_p75_4_3y, eigen_quad_p75_4_4y, eigen_quad_p75_4_5y, eigen_quad_p75_5_2y, eigen_quad_p75_5_3y, eigen_quad_p75_5_4y, eigen_quad_p75_5_5y, eigen_quad_p75_6_2y, eigen_quad_p75_6_3y, eigen_quad_p75_6_4y, eigen_quad_p75_6_5y, eigen_quad_p75_7_2y, eigen_quad_p75_7_3y, eigen_quad_p75_7_4y, eigen_quad_p75_7_5y, eigen_quad_p76_2_2y, eigen_quad_p76_2_3y, eigen_quad_p76_2_4y, eigen_quad_p76_2_5y, eigen_quad_p76_2_6y, eigen_quad_p76_3_2y, eigen_quad_p76_3_3y, eigen_quad_p76_3_4y, eigen_quad_p76_3_5y, eigen_quad_p76_3_6y, eigen_quad_p76_4_2y, eigen_quad_p76_4_3y, eigen_quad_p76_4_4y, eigen_quad_p76_4_5y, eigen_quad_p76_4_6y, eigen_quad_p76_5_2y, eigen_quad_p76_5_3y, eigen_quad_p76_5_4y, eigen_quad_p76_5_5y, eigen_quad_p76_5_6y, eigen_quad_p76_6_2y, eigen_quad_p76_6_3y, eigen_quad_p76_6_4y, eigen_quad_p76_6_5y, eigen_quad_p76_6_6y, eigen_quad_p76_7_2y, eigen_quad_p76_7_3y, eigen_quad_p76_7_4y, eigen_quad_p76_7_5y, eigen_quad_p76_7_6y, eigen_quad_p77_2_2y, eigen_quad_p77_2_3y, eigen_quad_p77_2_4y, eigen_quad_p77_2_5y, eigen_quad_p77_2_6y, eigen_quad_p77_2_7y, eigen_quad_p77_3_2y, eigen_quad_p77_3_3y, eigen_quad_p77_3_4y, eigen_quad_p77_3_5y, eigen_quad_p77_3_6y, eigen_quad_p77_3_7y, eigen_quad_p77_4_2y, eigen_quad_p77_4_3y, eigen_quad_p77_4_4y, eigen_quad_p77_4_5y, eigen_quad_p77_4_6y, eigen_quad_p77_4_7y, eigen_quad_p77_5_2y, eigen_quad_p77_5_3y, eigen_quad_p77_5_4y, eigen_quad_p77_5_5y, eigen_quad_p77_5_6y, eigen_quad_p77_5_7y, eigen_quad_p77_6_2y, eigen_quad_p77_6_3y, eigen_quad_p77_6_4y, eigen_quad_p77_6_5y, eigen_quad_p77_6_6y, eigen_quad_p77_6_7y, eigen_quad_p77_7_2y, eigen_quad_p77_7_3y, eigen_quad_p77_7_4y, eigen_quad_p77_7_5y, eigen_quad_p77_7_6y, eigen_quad_p77_7_7y, eigen_quad_p78_2_2y, eigen_quad_p78_2_3y, eigen_quad_p78_2_4y, eigen_quad_p78_2_5y, eigen_quad_p78_2_6y, eigen_quad_p78_2_7y, eigen_quad_p78_2_8y, eigen_quad_p78_3_2y, eigen_quad_p78_3_3y, eigen_quad_p78_3_4y, eigen_quad_p78_3_5y, eigen_quad_p78_3_6y, eigen_quad_p78_3_7y, eigen_quad_p78_3_8y, eigen_quad_p78_4_2y, eigen_quad_p78_4_3y, eigen_quad_p78_4_4y, eigen_quad_p78_4_5y, eigen_quad_p78_4_6y, eigen_quad_p78_4_7y, eigen_quad_p78_4_8y, eigen_quad_p78_5_2y, eigen_quad_p78_5_3y, eigen_quad_p78_5_4y, eigen_quad_p78_5_5y, eigen_quad_p78_5_6y, eigen_quad_p78_5_7y, eigen_quad_p78_5_8y, eigen_quad_p78_6_2y, eigen_quad_p78_6_3y, eigen_quad_p78_6_4y, eigen_quad_p78_6_5y, eigen_quad_p78_6_6y, eigen_quad_p78_6_7y, eigen_quad_p78_6_8y, eigen_quad_p78_7_2y, eigen_quad_p78_7_3y, eigen_quad_p78_7_4y, eigen_quad_p78_7_5y, eigen_quad_p78_7_6y, eigen_quad_p78_7_7y, eigen_quad_p78_7_8y, eigen_quad_p79_2_2y, eigen_quad_p79_2_3y, eigen_quad_p79_2_4y, eigen_quad_p79_2_5y, eigen_quad_p79_2_6y, eigen_quad_p79_2_7y, eigen_quad_p79_2_8y, eigen_quad_p79_2_9y, eigen_quad_p79_3_2y, eigen_quad_p79_3_3y, eigen_quad_p79_3_4y, eigen_quad_p79_3_5y, eigen_quad_p79_3_6y, eigen_quad_p79_3_7y, eigen_quad_p79_3_8y, eigen_quad_p79_3_9y, eigen_quad_p79_4_2y, eigen_quad_p79_4_3y, eigen_quad_p79_4_4y, eigen_quad_p79_4_5y, eigen_quad_p79_4_6y, eigen_quad_p79_4_7y, eigen_quad_p79_4_8y, eigen_quad_p79_4_9y, eigen_quad_p79_5_2y, eigen_quad_p79_5_3y, eigen_quad_p79_5_4y, eigen_quad_p79_5_5y, eigen_quad_p79_5_6y, eigen_quad_p79_5_7y, eigen_quad_p79_5_8y, eigen_quad_p79_5_9y, eigen_quad_p79_6_2y, eigen_quad_p79_6_3y, eigen_quad_p79_6_4y, eigen_quad_p79_6_5y, eigen_quad_p79_6_6y, eigen_quad_p79_6_7y, eigen_quad_p79_6_8y, eigen_quad_p79_6_9y, eigen_quad_p79_7_2y, eigen_quad_p79_7_3y, eigen_quad_p79_7_4y, eigen_quad_p79_7_5y, eigen_quad_p79_7_6y, eigen_quad_p79_7_7y, eigen_quad_p79_7_8y, eigen_quad_p79_7_9y, eigen_quad_p710_2_2y, eigen_quad_p710_2_3y, eigen_quad_p710_2_4y, eigen_quad_p710_2_5y, eigen_quad_p710_2_6y, eigen_quad_p710_2_7y, eigen_quad_p710_2_8y, eigen_quad_p710_2_9y, eigen_quad_p710_2_10y, eigen_quad_p710_3_2y, eigen_quad_p710_3_3y, eigen_quad_p710_3_4y, eigen_quad_p710_3_5y, eigen_quad_p710_3_6y, eigen_quad_p710_3_7y, eigen_quad_p710_3_8y, eigen_quad_p710_3_9y, eigen_quad_p710_3_10y, eigen_quad_p710_4_2y, eigen_quad_p710_4_3y, eigen_quad_p710_4_4y, eigen_quad_p710_4_5y, eigen_quad_p710_4_6y, eigen_quad_p710_4_7y, eigen_quad_p710_4_8y, eigen_quad_p710_4_9y, eigen_quad_p710_4_10y, eigen_quad_p710_5_2y, eigen_quad_p710_5_3y, eigen_quad_p710_5_4y, eigen_quad_p710_5_5y, eigen_quad_p710_5_6y, eigen_quad_p710_5_7y, eigen_quad_p710_5_8y, eigen_quad_p710_5_9y, eigen_quad_p710_5_10y, eigen_quad_p710_6_2y, eigen_quad_p710_6_3y, eigen_quad_p710_6_4y, eigen_quad_p710_6_5y, eigen_quad_p710_6_6y, eigen_quad_p710_6_7y, eigen_quad_p710_6_8y, eigen_quad_p710_6_9y, eigen_quad_p710_6_10y, eigen_quad_p710_7_2y, eigen_quad_p710_7_3y, eigen_quad_p710_7_4y, eigen_quad_p710_7_5y, eigen_quad_p710_7_6y, eigen_quad_p710_7_7y, eigen_quad_p710_7_8y, eigen_quad_p710_7_9y, eigen_quad_p710_7_10y,

  eigen_quad_p82_2_2y, eigen_quad_p82_3_2y, eigen_quad_p82_4_2y, eigen_quad_p82_5_2y, eigen_quad_p82_6_2y, eigen_quad_p82_7_2y, eigen_quad_p82_8_2y, eigen_quad_p83_2_2y, eigen_quad_p83_2_3y, eigen_quad_p83_3_2y, eigen_quad_p83_3_3y, eigen_quad_p83_4_2y, eigen_quad_p83_4_3y, eigen_quad_p83_5_2y, eigen_quad_p83_5_3y, eigen_quad_p83_6_2y, eigen_quad_p83_6_3y, eigen_quad_p83_7_2y, eigen_quad_p83_7_3y, eigen_quad_p83_8_2y, eigen_quad_p83_8_3y, eigen_quad_p84_2_2y, eigen_quad_p84_2_3y, eigen_quad_p84_2_4y, eigen_quad_p84_3_2y, eigen_quad_p84_3_3y, eigen_quad_p84_3_4y, eigen_quad_p84_4_2y, eigen_quad_p84_4_3y, eigen_quad_p84_4_4y, eigen_quad_p84_5_2y, eigen_quad_p84_5_3y, eigen_quad_p84_5_4y, eigen_quad_p84_6_2y, eigen_quad_p84_6_3y, eigen_quad_p84_6_4y, eigen_quad_p84_7_2y, eigen_quad_p84_7_3y, eigen_quad_p84_7_4y, eigen_quad_p84_8_2y, eigen_quad_p84_8_3y, eigen_quad_p84_8_4y, eigen_quad_p85_2_2y, eigen_quad_p85_2_3y, eigen_quad_p85_2_4y, eigen_quad_p85_2_5y, eigen_quad_p85_3_2y, eigen_quad_p85_3_3y, eigen_quad_p85_3_4y, eigen_quad_p85_3_5y, eigen_quad_p85_4_2y, eigen_quad_p85_4_3y, eigen_quad_p85_4_4y, eigen_quad_p85_4_5y, eigen_quad_p85_5_2y, eigen_quad_p85_5_3y, eigen_quad_p85_5_4y, eigen_quad_p85_5_5y, eigen_quad_p85_6_2y, eigen_quad_p85_6_3y, eigen_quad_p85_6_4y, eigen_quad_p85_6_5y, eigen_quad_p85_7_2y, eigen_quad_p85_7_3y, eigen_quad_p85_7_4y, eigen_quad_p85_7_5y, eigen_quad_p85_8_2y, eigen_quad_p85_8_3y, eigen_quad_p85_8_4y, eigen_quad_p85_8_5y, eigen_quad_p86_2_2y, eigen_quad_p86_2_3y, eigen_quad_p86_2_4y, eigen_quad_p86_2_5y, eigen_quad_p86_2_6y, eigen_quad_p86_3_2y, eigen_quad_p86_3_3y, eigen_quad_p86_3_4y, eigen_quad_p86_3_5y, eigen_quad_p86_3_6y, eigen_quad_p86_4_2y, eigen_quad_p86_4_3y, eigen_quad_p86_4_4y, eigen_quad_p86_4_5y, eigen_quad_p86_4_6y, eigen_quad_p86_5_2y, eigen_quad_p86_5_3y, eigen_quad_p86_5_4y, eigen_quad_p86_5_5y, eigen_quad_p86_5_6y, eigen_quad_p86_6_2y, eigen_quad_p86_6_3y, eigen_quad_p86_6_4y, eigen_quad_p86_6_5y, eigen_quad_p86_6_6y, eigen_quad_p86_7_2y, eigen_quad_p86_7_3y, eigen_quad_p86_7_4y, eigen_quad_p86_7_5y, eigen_quad_p86_7_6y, eigen_quad_p86_8_2y, eigen_quad_p86_8_3y, eigen_quad_p86_8_4y, eigen_quad_p86_8_5y, eigen_quad_p86_8_6y, eigen_quad_p87_2_2y, eigen_quad_p87_2_3y, eigen_quad_p87_2_4y, eigen_quad_p87_2_5y, eigen_quad_p87_2_6y, eigen_quad_p87_2_7y, eigen_quad_p87_3_2y, eigen_quad_p87_3_3y, eigen_quad_p87_3_4y, eigen_quad_p87_3_5y, eigen_quad_p87_3_6y, eigen_quad_p87_3_7y, eigen_quad_p87_4_2y, eigen_quad_p87_4_3y, eigen_quad_p87_4_4y, eigen_quad_p87_4_5y, eigen_quad_p87_4_6y, eigen_quad_p87_4_7y, eigen_quad_p87_5_2y, eigen_quad_p87_5_3y, eigen_quad_p87_5_4y, eigen_quad_p87_5_5y, eigen_quad_p87_5_6y, eigen_quad_p87_5_7y, eigen_quad_p87_6_2y, eigen_quad_p87_6_3y, eigen_quad_p87_6_4y, eigen_quad_p87_6_5y, eigen_quad_p87_6_6y, eigen_quad_p87_6_7y, eigen_quad_p87_7_2y, eigen_quad_p87_7_3y, eigen_quad_p87_7_4y, eigen_quad_p87_7_5y, eigen_quad_p87_7_6y, eigen_quad_p87_7_7y, eigen_quad_p87_8_2y, eigen_quad_p87_8_3y, eigen_quad_p87_8_4y, eigen_quad_p87_8_5y, eigen_quad_p87_8_6y, eigen_quad_p87_8_7y, eigen_quad_p88_2_2y, eigen_quad_p88_2_3y, eigen_quad_p88_2_4y, eigen_quad_p88_2_5y, eigen_quad_p88_2_6y, eigen_quad_p88_2_7y, eigen_quad_p88_2_8y, eigen_quad_p88_3_2y, eigen_quad_p88_3_3y, eigen_quad_p88_3_4y, eigen_quad_p88_3_5y, eigen_quad_p88_3_6y, eigen_quad_p88_3_7y, eigen_quad_p88_3_8y, eigen_quad_p88_4_2y, eigen_quad_p88_4_3y, eigen_quad_p88_4_4y, eigen_quad_p88_4_5y, eigen_quad_p88_4_6y, eigen_quad_p88_4_7y, eigen_quad_p88_4_8y, eigen_quad_p88_5_2y, eigen_quad_p88_5_3y, eigen_quad_p88_5_4y, eigen_quad_p88_5_5y, eigen_quad_p88_5_6y, eigen_quad_p88_5_7y, eigen_quad_p88_5_8y, eigen_quad_p88_6_2y, eigen_quad_p88_6_3y, eigen_quad_p88_6_4y, eigen_quad_p88_6_5y, eigen_quad_p88_6_6y, eigen_quad_p88_6_7y, eigen_quad_p88_6_8y, eigen_quad_p88_7_2y, eigen_quad_p88_7_3y, eigen_quad_p88_7_4y, eigen_quad_p88_7_5y, eigen_quad_p88_7_6y, eigen_quad_p88_7_7y, eigen_quad_p88_7_8y, eigen_quad_p88_8_2y, eigen_quad_p88_8_3y, eigen_quad_p88_8_4y, eigen_quad_p88_8_5y, eigen_quad_p88_8_6y, eigen_quad_p88_8_7y, eigen_quad_p88_8_8y, eigen_quad_p89_2_2y, eigen_quad_p89_2_3y, eigen_quad_p89_2_4y, eigen_quad_p89_2_5y, eigen_quad_p89_2_6y, eigen_quad_p89_2_7y, eigen_quad_p89_2_8y, eigen_quad_p89_2_9y, eigen_quad_p89_3_2y, eigen_quad_p89_3_3y, eigen_quad_p89_3_4y, eigen_quad_p89_3_5y, eigen_quad_p89_3_6y, eigen_quad_p89_3_7y, eigen_quad_p89_3_8y, eigen_quad_p89_3_9y, eigen_quad_p89_4_2y, eigen_quad_p89_4_3y, eigen_quad_p89_4_4y, eigen_quad_p89_4_5y, eigen_quad_p89_4_6y, eigen_quad_p89_4_7y, eigen_quad_p89_4_8y, eigen_quad_p89_4_9y, eigen_quad_p89_5_2y, eigen_quad_p89_5_3y, eigen_quad_p89_5_4y, eigen_quad_p89_5_5y, eigen_quad_p89_5_6y, eigen_quad_p89_5_7y, eigen_quad_p89_5_8y, eigen_quad_p89_5_9y, eigen_quad_p89_6_2y, eigen_quad_p89_6_3y, eigen_quad_p89_6_4y, eigen_quad_p89_6_5y, eigen_quad_p89_6_6y, eigen_quad_p89_6_7y, eigen_quad_p89_6_8y, eigen_quad_p89_6_9y, eigen_quad_p89_7_2y, eigen_quad_p89_7_3y, eigen_quad_p89_7_4y, eigen_quad_p89_7_5y, eigen_quad_p89_7_6y, eigen_quad_p89_7_7y, eigen_quad_p89_7_8y, eigen_quad_p89_7_9y, eigen_quad_p89_8_2y, eigen_quad_p89_8_3y, eigen_quad_p89_8_4y, eigen_quad_p89_8_5y, eigen_quad_p89_8_6y, eigen_quad_p89_8_7y, eigen_quad_p89_8_8y, eigen_quad_p89_8_9y, eigen_quad_p810_2_2y, eigen_quad_p810_2_3y, eigen_quad_p810_2_4y, eigen_quad_p810_2_5y, eigen_quad_p810_2_6y, eigen_quad_p810_2_7y, eigen_quad_p810_2_8y, eigen_quad_p810_2_9y, eigen_quad_p810_2_10y, eigen_quad_p810_3_2y, eigen_quad_p810_3_3y, eigen_quad_p810_3_4y, eigen_quad_p810_3_5y, eigen_quad_p810_3_6y, eigen_quad_p810_3_7y, eigen_quad_p810_3_8y, eigen_quad_p810_3_9y, eigen_quad_p810_3_10y, eigen_quad_p810_4_2y, eigen_quad_p810_4_3y, eigen_quad_p810_4_4y, eigen_quad_p810_4_5y, eigen_quad_p810_4_6y, eigen_quad_p810_4_7y, eigen_quad_p810_4_8y, eigen_quad_p810_4_9y, eigen_quad_p810_4_10y, eigen_quad_p810_5_2y, eigen_quad_p810_5_3y, eigen_quad_p810_5_4y, eigen_quad_p810_5_5y, eigen_quad_p810_5_6y, eigen_quad_p810_5_7y, eigen_quad_p810_5_8y, eigen_quad_p810_5_9y, eigen_quad_p810_5_10y, eigen_quad_p810_6_2y, eigen_quad_p810_6_3y, eigen_quad_p810_6_4y, eigen_quad_p810_6_5y, eigen_quad_p810_6_6y, eigen_quad_p810_6_7y, eigen_quad_p810_6_8y, eigen_quad_p810_6_9y, eigen_quad_p810_6_10y, eigen_quad_p810_7_2y, eigen_quad_p810_7_3y, eigen_quad_p810_7_4y, eigen_quad_p810_7_5y, eigen_quad_p810_7_6y, eigen_quad_p810_7_7y, eigen_quad_p810_7_8y, eigen_quad_p810_7_9y, eigen_quad_p810_7_10y, eigen_quad_p810_8_2y, eigen_quad_p810_8_3y, eigen_quad_p810_8_4y, eigen_quad_p810_8_5y, eigen_quad_p810_8_6y, eigen_quad_p810_8_7y, eigen_quad_p810_8_8y, eigen_quad_p810_8_9y, eigen_quad_p810_8_10y,

  eigen_quad_p92_2_2y, eigen_quad_p92_3_2y, eigen_quad_p92_4_2y, eigen_quad_p92_5_2y, eigen_quad_p92_6_2y, eigen_quad_p92_7_2y, eigen_quad_p92_8_2y, eigen_quad_p92_9_2y, eigen_quad_p93_2_2y, eigen_quad_p93_2_3y, eigen_quad_p93_3_2y, eigen_quad_p93_3_3y, eigen_quad_p93_4_2y, eigen_quad_p93_4_3y, eigen_quad_p93_5_2y, eigen_quad_p93_5_3y, eigen_quad_p93_6_2y, eigen_quad_p93_6_3y, eigen_quad_p93_7_2y, eigen_quad_p93_7_3y, eigen_quad_p93_8_2y, eigen_quad_p93_8_3y, eigen_quad_p93_9_2y, eigen_quad_p93_9_3y, eigen_quad_p94_2_2y, eigen_quad_p94_2_3y, eigen_quad_p94_2_4y, eigen_quad_p94_3_2y, eigen_quad_p94_3_3y, eigen_quad_p94_3_4y, eigen_quad_p94_4_2y, eigen_quad_p94_4_3y, eigen_quad_p94_4_4y, eigen_quad_p94_5_2y, eigen_quad_p94_5_3y, eigen_quad_p94_5_4y, eigen_quad_p94_6_2y, eigen_quad_p94_6_3y, eigen_quad_p94_6_4y, eigen_quad_p94_7_2y, eigen_quad_p94_7_3y, eigen_quad_p94_7_4y, eigen_quad_p94_8_2y, eigen_quad_p94_8_3y, eigen_quad_p94_8_4y, eigen_quad_p94_9_2y, eigen_quad_p94_9_3y, eigen_quad_p94_9_4y, eigen_quad_p95_2_2y, eigen_quad_p95_2_3y, eigen_quad_p95_2_4y, eigen_quad_p95_2_5y, eigen_quad_p95_3_2y, eigen_quad_p95_3_3y, eigen_quad_p95_3_4y, eigen_quad_p95_3_5y, eigen_quad_p95_4_2y, eigen_quad_p95_4_3y, eigen_quad_p95_4_4y, eigen_quad_p95_4_5y, eigen_quad_p95_5_2y, eigen_quad_p95_5_3y, eigen_quad_p95_5_4y, eigen_quad_p95_5_5y, eigen_quad_p95_6_2y, eigen_quad_p95_6_3y, eigen_quad_p95_6_4y, eigen_quad_p95_6_5y, eigen_quad_p95_7_2y, eigen_quad_p95_7_3y, eigen_quad_p95_7_4y, eigen_quad_p95_7_5y, eigen_quad_p95_8_2y, eigen_quad_p95_8_3y, eigen_quad_p95_8_4y, eigen_quad_p95_8_5y, eigen_quad_p95_9_2y, eigen_quad_p95_9_3y, eigen_quad_p95_9_4y, eigen_quad_p95_9_5y, eigen_quad_p96_2_2y, eigen_quad_p96_2_3y, eigen_quad_p96_2_4y, eigen_quad_p96_2_5y, eigen_quad_p96_2_6y, eigen_quad_p96_3_2y, eigen_quad_p96_3_3y, eigen_quad_p96_3_4y, eigen_quad_p96_3_5y, eigen_quad_p96_3_6y, eigen_quad_p96_4_2y, eigen_quad_p96_4_3y, eigen_quad_p96_4_4y, eigen_quad_p96_4_5y, eigen_quad_p96_4_6y, eigen_quad_p96_5_2y, eigen_quad_p96_5_3y, eigen_quad_p96_5_4y, eigen_quad_p96_5_5y, eigen_quad_p96_5_6y, eigen_quad_p96_6_2y, eigen_quad_p96_6_3y, eigen_quad_p96_6_4y, eigen_quad_p96_6_5y, eigen_quad_p96_6_6y, eigen_quad_p96_7_2y, eigen_quad_p96_7_3y, eigen_quad_p96_7_4y, eigen_quad_p96_7_5y, eigen_quad_p96_7_6y, eigen_quad_p96_8_2y, eigen_quad_p96_8_3y, eigen_quad_p96_8_4y, eigen_quad_p96_8_5y, eigen_quad_p96_8_6y, eigen_quad_p96_9_2y, eigen_quad_p96_9_3y, eigen_quad_p96_9_4y, eigen_quad_p96_9_5y, eigen_quad_p96_9_6y, eigen_quad_p97_2_2y, eigen_quad_p97_2_3y, eigen_quad_p97_2_4y, eigen_quad_p97_2_5y, eigen_quad_p97_2_6y, eigen_quad_p97_2_7y, eigen_quad_p97_3_2y, eigen_quad_p97_3_3y, eigen_quad_p97_3_4y, eigen_quad_p97_3_5y, eigen_quad_p97_3_6y, eigen_quad_p97_3_7y, eigen_quad_p97_4_2y, eigen_quad_p97_4_3y, eigen_quad_p97_4_4y, eigen_quad_p97_4_5y, eigen_quad_p97_4_6y, eigen_quad_p97_4_7y, eigen_quad_p97_5_2y, eigen_quad_p97_5_3y, eigen_quad_p97_5_4y, eigen_quad_p97_5_5y, eigen_quad_p97_5_6y, eigen_quad_p97_5_7y, eigen_quad_p97_6_2y, eigen_quad_p97_6_3y, eigen_quad_p97_6_4y, eigen_quad_p97_6_5y, eigen_quad_p97_6_6y, eigen_quad_p97_6_7y, eigen_quad_p97_7_2y, eigen_quad_p97_7_3y, eigen_quad_p97_7_4y, eigen_quad_p97_7_5y, eigen_quad_p97_7_6y, eigen_quad_p97_7_7y, eigen_quad_p97_8_2y, eigen_quad_p97_8_3y, eigen_quad_p97_8_4y, eigen_quad_p97_8_5y, eigen_quad_p97_8_6y, eigen_quad_p97_8_7y, eigen_quad_p97_9_2y, eigen_quad_p97_9_3y, eigen_quad_p97_9_4y, eigen_quad_p97_9_5y, eigen_quad_p97_9_6y, eigen_quad_p97_9_7y, eigen_quad_p98_2_2y, eigen_quad_p98_2_3y, eigen_quad_p98_2_4y, eigen_quad_p98_2_5y, eigen_quad_p98_2_6y, eigen_quad_p98_2_7y, eigen_quad_p98_2_8y, eigen_quad_p98_3_2y, eigen_quad_p98_3_3y, eigen_quad_p98_3_4y, eigen_quad_p98_3_5y, eigen_quad_p98_3_6y, eigen_quad_p98_3_7y, eigen_quad_p98_3_8y, eigen_quad_p98_4_2y, eigen_quad_p98_4_3y, eigen_quad_p98_4_4y, eigen_quad_p98_4_5y, eigen_quad_p98_4_6y, eigen_quad_p98_4_7y, eigen_quad_p98_4_8y, eigen_quad_p98_5_2y, eigen_quad_p98_5_3y, eigen_quad_p98_5_4y, eigen_quad_p98_5_5y, eigen_quad_p98_5_6y, eigen_quad_p98_5_7y, eigen_quad_p98_5_8y, eigen_quad_p98_6_2y, eigen_quad_p98_6_3y, eigen_quad_p98_6_4y, eigen_quad_p98_6_5y, eigen_quad_p98_6_6y, eigen_quad_p98_6_7y, eigen_quad_p98_6_8y, eigen_quad_p98_7_2y, eigen_quad_p98_7_3y, eigen_quad_p98_7_4y, eigen_quad_p98_7_5y, eigen_quad_p98_7_6y, eigen_quad_p98_7_7y, eigen_quad_p98_7_8y, eigen_quad_p98_8_2y, eigen_quad_p98_8_3y, eigen_quad_p98_8_4y, eigen_quad_p98_8_5y, eigen_quad_p98_8_6y, eigen_quad_p98_8_7y, eigen_quad_p98_8_8y, eigen_quad_p98_9_2y, eigen_quad_p98_9_3y, eigen_quad_p98_9_4y, eigen_quad_p98_9_5y, eigen_quad_p98_9_6y, eigen_quad_p98_9_7y, eigen_quad_p98_9_8y, eigen_quad_p99_2_2y, eigen_quad_p99_2_3y, eigen_quad_p99_2_4y, eigen_quad_p99_2_5y, eigen_quad_p99_2_6y, eigen_quad_p99_2_7y, eigen_quad_p99_2_8y, eigen_quad_p99_2_9y, eigen_quad_p99_3_2y, eigen_quad_p99_3_3y, eigen_quad_p99_3_4y, eigen_quad_p99_3_5y, eigen_quad_p99_3_6y, eigen_quad_p99_3_7y, eigen_quad_p99_3_8y, eigen_quad_p99_3_9y, eigen_quad_p99_4_2y, eigen_quad_p99_4_3y, eigen_quad_p99_4_4y, eigen_quad_p99_4_5y, eigen_quad_p99_4_6y, eigen_quad_p99_4_7y, eigen_quad_p99_4_8y, eigen_quad_p99_4_9y, eigen_quad_p99_5_2y, eigen_quad_p99_5_3y, eigen_quad_p99_5_4y, eigen_quad_p99_5_5y, eigen_quad_p99_5_6y, eigen_quad_p99_5_7y, eigen_quad_p99_5_8y, eigen_quad_p99_5_9y, eigen_quad_p99_6_2y, eigen_quad_p99_6_3y, eigen_quad_p99_6_4y, eigen_quad_p99_6_5y, eigen_quad_p99_6_6y, eigen_quad_p99_6_7y, eigen_quad_p99_6_8y, eigen_quad_p99_6_9y, eigen_quad_p99_7_2y, eigen_quad_p99_7_3y, eigen_quad_p99_7_4y, eigen_quad_p99_7_5y, eigen_quad_p99_7_6y, eigen_quad_p99_7_7y, eigen_quad_p99_7_8y, eigen_quad_p99_7_9y, eigen_quad_p99_8_2y, eigen_quad_p99_8_3y, eigen_quad_p99_8_4y, eigen_quad_p99_8_5y, eigen_quad_p99_8_6y, eigen_quad_p99_8_7y, eigen_quad_p99_8_8y, eigen_quad_p99_8_9y, eigen_quad_p99_9_2y, eigen_quad_p99_9_3y, eigen_quad_p99_9_4y, eigen_quad_p99_9_5y, eigen_quad_p99_9_6y, eigen_quad_p99_9_7y, eigen_quad_p99_9_8y, eigen_quad_p99_9_9y, eigen_quad_p910_2_2y, eigen_quad_p910_2_3y, eigen_quad_p910_2_4y, eigen_quad_p910_2_5y, eigen_quad_p910_2_6y, eigen_quad_p910_2_7y, eigen_quad_p910_2_8y, eigen_quad_p910_2_9y, eigen_quad_p910_2_10y, eigen_quad_p910_3_2y, eigen_quad_p910_3_3y, eigen_quad_p910_3_4y, eigen_quad_p910_3_5y, eigen_quad_p910_3_6y, eigen_quad_p910_3_7y, eigen_quad_p910_3_8y, eigen_quad_p910_3_9y, eigen_quad_p910_3_10y, eigen_quad_p910_4_2y, eigen_quad_p910_4_3y, eigen_quad_p910_4_4y, eigen_quad_p910_4_5y, eigen_quad_p910_4_6y, eigen_quad_p910_4_7y, eigen_quad_p910_4_8y, eigen_quad_p910_4_9y, eigen_quad_p910_4_10y, eigen_quad_p910_5_2y, eigen_quad_p910_5_3y, eigen_quad_p910_5_4y, eigen_quad_p910_5_5y, eigen_quad_p910_5_6y, eigen_quad_p910_5_7y, eigen_quad_p910_5_8y, eigen_quad_p910_5_9y, eigen_quad_p910_5_10y, eigen_quad_p910_6_2y, eigen_quad_p910_6_3y, eigen_quad_p910_6_4y, eigen_quad_p910_6_5y, eigen_quad_p910_6_6y, eigen_quad_p910_6_7y, eigen_quad_p910_6_8y, eigen_quad_p910_6_9y, eigen_quad_p910_6_10y, eigen_quad_p910_7_2y, eigen_quad_p910_7_3y, eigen_quad_p910_7_4y, eigen_quad_p910_7_5y, eigen_quad_p910_7_6y, eigen_quad_p910_7_7y, eigen_quad_p910_7_8y, eigen_quad_p910_7_9y, eigen_quad_p910_7_10y, eigen_quad_p910_8_2y, eigen_quad_p910_8_3y, eigen_quad_p910_8_4y, eigen_quad_p910_8_5y, eigen_quad_p910_8_6y, eigen_quad_p910_8_7y, eigen_quad_p910_8_8y, eigen_quad_p910_8_9y, eigen_quad_p910_8_10y, eigen_quad_p910_9_2y, eigen_quad_p910_9_3y, eigen_quad_p910_9_4y, eigen_quad_p910_9_5y, eigen_quad_p910_9_6y, eigen_quad_p910_9_7y, eigen_quad_p910_9_8y, eigen_quad_p910_9_9y, eigen_quad_p910_9_10y,

  eigen_quad_p102_2_2y, eigen_quad_p102_3_2y, eigen_quad_p102_4_2y, eigen_quad_p102_5_2y, eigen_quad_p102_6_2y, eigen_quad_p102_7_2y, eigen_quad_p102_8_2y, eigen_quad_p102_9_2y, eigen_quad_p102_10_2y, eigen_quad_p103_2_2y, eigen_quad_p103_2_3y, eigen_quad_p103_3_2y, eigen_quad_p103_3_3y, eigen_quad_p103_4_2y, eigen_quad_p103_4_3y, eigen_quad_p103_5_2y, eigen_quad_p103_5_3y, eigen_quad_p103_6_2y, eigen_quad_p103_6_3y, eigen_quad_p103_7_2y, eigen_quad_p103_7_3y, eigen_quad_p103_8_2y, eigen_quad_p103_8_3y, eigen_quad_p103_9_2y, eigen_quad_p103_9_3y, eigen_quad_p103_10_2y, eigen_quad_p103_10_3y, eigen_quad_p104_2_2y, eigen_quad_p104_2_3y, eigen_quad_p104_2_4y, eigen_quad_p104_3_2y, eigen_quad_p104_3_3y, eigen_quad_p104_3_4y, eigen_quad_p104_4_2y, eigen_quad_p104_4_3y, eigen_quad_p104_4_4y, eigen_quad_p104_5_2y, eigen_quad_p104_5_3y, eigen_quad_p104_5_4y, eigen_quad_p104_6_2y, eigen_quad_p104_6_3y, eigen_quad_p104_6_4y, eigen_quad_p104_7_2y, eigen_quad_p104_7_3y, eigen_quad_p104_7_4y, eigen_quad_p104_8_2y, eigen_quad_p104_8_3y, eigen_quad_p104_8_4y, eigen_quad_p104_9_2y, eigen_quad_p104_9_3y, eigen_quad_p104_9_4y, eigen_quad_p104_10_2y, eigen_quad_p104_10_3y, eigen_quad_p104_10_4y, eigen_quad_p105_2_2y, eigen_quad_p105_2_3y, eigen_quad_p105_2_4y, eigen_quad_p105_2_5y, eigen_quad_p105_3_2y, eigen_quad_p105_3_3y, eigen_quad_p105_3_4y, eigen_quad_p105_3_5y, eigen_quad_p105_4_2y, eigen_quad_p105_4_3y, eigen_quad_p105_4_4y, eigen_quad_p105_4_5y, eigen_quad_p105_5_2y, eigen_quad_p105_5_3y, eigen_quad_p105_5_4y, eigen_quad_p105_5_5y, eigen_quad_p105_6_2y, eigen_quad_p105_6_3y, eigen_quad_p105_6_4y, eigen_quad_p105_6_5y, eigen_quad_p105_7_2y, eigen_quad_p105_7_3y, eigen_quad_p105_7_4y, eigen_quad_p105_7_5y, eigen_quad_p105_8_2y, eigen_quad_p105_8_3y, eigen_quad_p105_8_4y, eigen_quad_p105_8_5y, eigen_quad_p105_9_2y, eigen_quad_p105_9_3y, eigen_quad_p105_9_4y, eigen_quad_p105_9_5y, eigen_quad_p105_10_2y, eigen_quad_p105_10_3y, eigen_quad_p105_10_4y, eigen_quad_p105_10_5y, eigen_quad_p106_2_2y, eigen_quad_p106_2_3y, eigen_quad_p106_2_4y, eigen_quad_p106_2_5y, eigen_quad_p106_2_6y, eigen_quad_p106_3_2y, eigen_quad_p106_3_3y, eigen_quad_p106_3_4y, eigen_quad_p106_3_5y, eigen_quad_p106_3_6y, eigen_quad_p106_4_2y, eigen_quad_p106_4_3y, eigen_quad_p106_4_4y, eigen_quad_p106_4_5y, eigen_quad_p106_4_6y, eigen_quad_p106_5_2y, eigen_quad_p106_5_3y, eigen_quad_p106_5_4y, eigen_quad_p106_5_5y, eigen_quad_p106_5_6y, eigen_quad_p106_6_2y, eigen_quad_p106_6_3y, eigen_quad_p106_6_4y, eigen_quad_p106_6_5y, eigen_quad_p106_6_6y, eigen_quad_p106_7_2y, eigen_quad_p106_7_3y, eigen_quad_p106_7_4y, eigen_quad_p106_7_5y, eigen_quad_p106_7_6y, eigen_quad_p106_8_2y, eigen_quad_p106_8_3y, eigen_quad_p106_8_4y, eigen_quad_p106_8_5y, eigen_quad_p106_8_6y, eigen_quad_p106_9_2y, eigen_quad_p106_9_3y, eigen_quad_p106_9_4y, eigen_quad_p106_9_5y, eigen_quad_p106_9_6y, eigen_quad_p106_10_2y, eigen_quad_p106_10_3y, eigen_quad_p106_10_4y, eigen_quad_p106_10_5y, eigen_quad_p106_10_6y, eigen_quad_p107_2_2y, eigen_quad_p107_2_3y, eigen_quad_p107_2_4y, eigen_quad_p107_2_5y, eigen_quad_p107_2_6y, eigen_quad_p107_2_7y, eigen_quad_p107_3_2y, eigen_quad_p107_3_3y, eigen_quad_p107_3_4y, eigen_quad_p107_3_5y, eigen_quad_p107_3_6y, eigen_quad_p107_3_7y, eigen_quad_p107_4_2y, eigen_quad_p107_4_3y, eigen_quad_p107_4_4y, eigen_quad_p107_4_5y, eigen_quad_p107_4_6y, eigen_quad_p107_4_7y, eigen_quad_p107_5_2y, eigen_quad_p107_5_3y, eigen_quad_p107_5_4y, eigen_quad_p107_5_5y, eigen_quad_p107_5_6y, eigen_quad_p107_5_7y, eigen_quad_p107_6_2y, eigen_quad_p107_6_3y, eigen_quad_p107_6_4y, eigen_quad_p107_6_5y, eigen_quad_p107_6_6y, eigen_quad_p107_6_7y, eigen_quad_p107_7_2y, eigen_quad_p107_7_3y, eigen_quad_p107_7_4y, eigen_quad_p107_7_5y, eigen_quad_p107_7_6y, eigen_quad_p107_7_7y, eigen_quad_p107_8_2y, eigen_quad_p107_8_3y, eigen_quad_p107_8_4y, eigen_quad_p107_8_5y, eigen_quad_p107_8_6y, eigen_quad_p107_8_7y, eigen_quad_p107_9_2y, eigen_quad_p107_9_3y, eigen_quad_p107_9_4y, eigen_quad_p107_9_5y, eigen_quad_p107_9_6y, eigen_quad_p107_9_7y, eigen_quad_p107_10_2y, eigen_quad_p107_10_3y, eigen_quad_p107_10_4y, eigen_quad_p107_10_5y, eigen_quad_p107_10_6y, eigen_quad_p107_10_7y, eigen_quad_p108_2_2y, eigen_quad_p108_2_3y, eigen_quad_p108_2_4y, eigen_quad_p108_2_5y, eigen_quad_p108_2_6y, eigen_quad_p108_2_7y, eigen_quad_p108_2_8y, eigen_quad_p108_3_2y, eigen_quad_p108_3_3y, eigen_quad_p108_3_4y, eigen_quad_p108_3_5y, eigen_quad_p108_3_6y, eigen_quad_p108_3_7y, eigen_quad_p108_3_8y, eigen_quad_p108_4_2y, eigen_quad_p108_4_3y, eigen_quad_p108_4_4y, eigen_quad_p108_4_5y, eigen_quad_p108_4_6y, eigen_quad_p108_4_7y, eigen_quad_p108_4_8y, eigen_quad_p108_5_2y, eigen_quad_p108_5_3y, eigen_quad_p108_5_4y, eigen_quad_p108_5_5y, eigen_quad_p108_5_6y, eigen_quad_p108_5_7y, eigen_quad_p108_5_8y, eigen_quad_p108_6_2y, eigen_quad_p108_6_3y, eigen_quad_p108_6_4y, eigen_quad_p108_6_5y, eigen_quad_p108_6_6y, eigen_quad_p108_6_7y, eigen_quad_p108_6_8y, eigen_quad_p108_7_2y, eigen_quad_p108_7_3y, eigen_quad_p108_7_4y, eigen_quad_p108_7_5y, eigen_quad_p108_7_6y, eigen_quad_p108_7_7y, eigen_quad_p108_7_8y, eigen_quad_p108_8_2y, eigen_quad_p108_8_3y, eigen_quad_p108_8_4y, eigen_quad_p108_8_5y, eigen_quad_p108_8_6y, eigen_quad_p108_8_7y, eigen_quad_p108_8_8y, eigen_quad_p108_9_2y, eigen_quad_p108_9_3y, eigen_quad_p108_9_4y, eigen_quad_p108_9_5y, eigen_quad_p108_9_6y, eigen_quad_p108_9_7y, eigen_quad_p108_9_8y, eigen_quad_p108_10_2y, eigen_quad_p108_10_3y, eigen_quad_p108_10_4y, eigen_quad_p108_10_5y, eigen_quad_p108_10_6y, eigen_quad_p108_10_7y, eigen_quad_p108_10_8y, eigen_quad_p109_2_2y, eigen_quad_p109_2_3y, eigen_quad_p109_2_4y, eigen_quad_p109_2_5y, eigen_quad_p109_2_6y, eigen_quad_p109_2_7y, eigen_quad_p109_2_8y, eigen_quad_p109_2_9y, eigen_quad_p109_3_2y, eigen_quad_p109_3_3y, eigen_quad_p109_3_4y, eigen_quad_p109_3_5y, eigen_quad_p109_3_6y, eigen_quad_p109_3_7y, eigen_quad_p109_3_8y, eigen_quad_p109_3_9y, eigen_quad_p109_4_2y, eigen_quad_p109_4_3y, eigen_quad_p109_4_4y, eigen_quad_p109_4_5y, eigen_quad_p109_4_6y, eigen_quad_p109_4_7y, eigen_quad_p109_4_8y, eigen_quad_p109_4_9y, eigen_quad_p109_5_2y, eigen_quad_p109_5_3y, eigen_quad_p109_5_4y, eigen_quad_p109_5_5y, eigen_quad_p109_5_6y, eigen_quad_p109_5_7y, eigen_quad_p109_5_8y, eigen_quad_p109_5_9y, eigen_quad_p109_6_2y, eigen_quad_p109_6_3y, eigen_quad_p109_6_4y, eigen_quad_p109_6_5y, eigen_quad_p109_6_6y, eigen_quad_p109_6_7y, eigen_quad_p109_6_8y, eigen_quad_p109_6_9y, eigen_quad_p109_7_2y, eigen_quad_p109_7_3y, eigen_quad_p109_7_4y, eigen_quad_p109_7_5y, eigen_quad_p109_7_6y, eigen_quad_p109_7_7y, eigen_quad_p109_7_8y, eigen_quad_p109_7_9y, eigen_quad_p109_8_2y, eigen_quad_p109_8_3y, eigen_quad_p109_8_4y, eigen_quad_p109_8_5y, eigen_quad_p109_8_6y, eigen_quad_p109_8_7y, eigen_quad_p109_8_8y, eigen_quad_p109_8_9y, eigen_quad_p109_9_2y, eigen_quad_p109_9_3y, eigen_quad_p109_9_4y, eigen_quad_p109_9_5y, eigen_quad_p109_9_6y, eigen_quad_p109_9_7y, eigen_quad_p109_9_8y, eigen_quad_p109_9_9y, eigen_quad_p109_10_2y, eigen_quad_p109_10_3y, eigen_quad_p109_10_4y, eigen_quad_p109_10_5y, eigen_quad_p109_10_6y, eigen_quad_p109_10_7y, eigen_quad_p109_10_8y, eigen_quad_p109_10_9y, eigen_quad_p1010_2_2y, eigen_quad_p1010_2_3y, eigen_quad_p1010_2_4y, eigen_quad_p1010_2_5y, eigen_quad_p1010_2_6y, eigen_quad_p1010_2_7y, eigen_quad_p1010_2_8y, eigen_quad_p1010_2_9y, eigen_quad_p1010_2_10y, eigen_quad_p1010_3_2y, eigen_quad_p1010_3_3y, eigen_quad_p1010_3_4y, eigen_quad_p1010_3_5y, eigen_quad_p1010_3_6y, eigen_quad_p1010_3_7y, eigen_quad_p1010_3_8y, eigen_quad_p1010_3_9y, eigen_quad_p1010_3_10y, eigen_quad_p1010_4_2y, eigen_quad_p1010_4_3y, eigen_quad_p1010_4_4y, eigen_quad_p1010_4_5y, eigen_quad_p1010_4_6y, eigen_quad_p1010_4_7y, eigen_quad_p1010_4_8y, eigen_quad_p1010_4_9y, eigen_quad_p1010_4_10y, eigen_quad_p1010_5_2y, eigen_quad_p1010_5_3y, eigen_quad_p1010_5_4y, eigen_quad_p1010_5_5y, eigen_quad_p1010_5_6y, eigen_quad_p1010_5_7y, eigen_quad_p1010_5_8y, eigen_quad_p1010_5_9y, eigen_quad_p1010_5_10y, eigen_quad_p1010_6_2y, eigen_quad_p1010_6_3y, eigen_quad_p1010_6_4y, eigen_quad_p1010_6_5y, eigen_quad_p1010_6_6y, eigen_quad_p1010_6_7y, eigen_quad_p1010_6_8y, eigen_quad_p1010_6_9y, eigen_quad_p1010_6_10y, eigen_quad_p1010_7_2y, eigen_quad_p1010_7_3y, eigen_quad_p1010_7_4y, eigen_quad_p1010_7_5y, eigen_quad_p1010_7_6y, eigen_quad_p1010_7_7y, eigen_quad_p1010_7_8y, eigen_quad_p1010_7_9y, eigen_quad_p1010_7_10y, eigen_quad_p1010_8_2y, eigen_quad_p1010_8_3y, eigen_quad_p1010_8_4y, eigen_quad_p1010_8_5y, eigen_quad_p1010_8_6y, eigen_quad_p1010_8_7y, eigen_quad_p1010_8_8y, eigen_quad_p1010_8_9y, eigen_quad_p1010_8_10y, eigen_quad_p1010_9_2y, eigen_quad_p1010_9_3y, eigen_quad_p1010_9_4y, eigen_quad_p1010_9_5y, eigen_quad_p1010_9_6y, eigen_quad_p1010_9_7y, eigen_quad_p1010_9_8y, eigen_quad_p1010_9_9y, eigen_quad_p1010_9_10y, eigen_quad_p1010_10_2y, eigen_quad_p1010_10_3y, eigen_quad_p1010_10_4y, eigen_quad_p1010_10_5y, eigen_quad_p1010_10_6y, eigen_quad_p1010_10_7y, eigen_quad_p1010_10_8y, eigen_quad_p1010_10_9y, eigen_quad_p1010_10_10y
};


  static int qb_2_2[] = { 56,};
  static int qb_2_3[] = { 57,58,};
  static int qb_2_4[] = { 59,60,61,};
  static int qb_2_5[] = { 62,63,64,65,};
  static int qb_2_6[] = { 66,67,68,69,70,};
  static int qb_2_7[] = { 71,72,73,74,75,76,};
  static int qb_2_8[] = { 77,78,79,80,81,82,83,};
  static int qb_2_9[] = { 84,85,86,87,88,89,90,91,};
  static int qb_2_10[] = { 92,93,94,95,96,97,98,99,100,};

  static int qb_3_2[] = { 101,102,};
  static int qb_3_3[] = { 103,104,105,106,};
  static int qb_3_4[] = { 107,108,109,110,111,112,};
  static int qb_3_5[] = { 113,114,115,116,117,118,119,120,};
  static int qb_3_6[] = { 121,122,123,124,125,126,127,128,129,130,};
  static int qb_3_7[] = { 131,132,133,134,135,136,137,138,139,140,141,142,};
  static int qb_3_8[] = { 143,144,145,146,147,148,149,150,151,152,153,154,155,156,};
  static int qb_3_9[] = { 157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,};
  static int qb_3_10[] = { 173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,};

  static int qb_4_2[] = { 191,192,193,};
  static int qb_4_3[] = { 194,195,196,197,198,199,};
  static int qb_4_4[] = { 200,201,202,203,204,205,206,207,208,};
  static int qb_4_5[] = { 209,210,211,212,213,214,215,216,217,218,219,220,};
  static int qb_4_6[] = { 221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,};
  static int qb_4_7[] = { 236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,};
  static int qb_4_8[] = { 254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,};
  static int qb_4_9[] = { 275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,};
  static int qb_4_10[] = { 299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,};

  static int qb_5_2[] = { 326,327,328,329,};
  static int qb_5_3[] = { 330,331,332,333,334,335,336,337,};
  static int qb_5_4[] = { 338,339,340,341,342,343,344,345,346,347,348,349,};
  static int qb_5_5[] = { 350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,};
  static int qb_5_6[] = { 366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,};
  static int qb_5_7[] = { 386,387,388,389,390,391,392,393,394,395,396,397,398,399,400,401,402,403,404,405,406,407,408,409,};
  static int qb_5_8[] = { 410,411,412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428,429,430,431,432,433,434,435,436,437,};
  static int qb_5_9[] = { 438,439,440,441,442,443,444,445,446,447,448,449,450,451,452,453,454,455,456,457,458,459,460,461,462,463,464,465,466,467,468,469,};
  static int qb_5_10[] = { 470,471,472,473,474,475,476,477,478,479,480,481,482,483,484,485,486,487,488,489,490,491,492,493,494,495,496,497,498,499,500,501,502,503,504,505,};

  static int qb_6_2[] = { 506,507,508,509,510,};
  static int qb_6_3[] = { 511,512,513,514,515,516,517,518,519,520,};
  static int qb_6_4[] = { 521,522,523,524,525,526,527,528,529,530,531,532,533,534,535,};
  static int qb_6_5[] = { 536,537,538,539,540,541,542,543,544,545,546,547,548,549,550,551,552,553,554,555,};
  static int qb_6_6[] = { 556,557,558,559,560,561,562,563,564,565,566,567,568,569,570,571,572,573,574,575,576,577,578,579,580,};
  static int qb_6_7[] = { 581,582,583,584,585,586,587,588,589,590,591,592,593,594,595,596,597,598,599,600,601,602,603,604,605,606,607,608,609,610,};
  static int qb_6_8[] = { 611,612,613,614,615,616,617,618,619,620,621,622,623,624,625,626,627,628,629,630,631,632,633,634,635,636,637,638,639,640,641,642,643,644,645,};
  static int qb_6_9[] = { 646,647,648,649,650,651,652,653,654,655,656,657,658,659,660,661,662,663,664,665,666,667,668,669,670,671,672,673,674,675,676,677,678,679,680,681,682,683,684,685,};
  static int qb_6_10[] = { 686,687,688,689,690,691,692,693,694,695,696,697,698,699,700,701,702,703,704,705,706,707,708,709,710,711,712,713,714,715,716,717,718,719,720,721,722,723,724,725,726,727,728,729,730,};

  static int qb_7_2[] = { 731,732,733,734,735,736,};
  static int qb_7_3[] = { 737,738,739,740,741,742,743,744,745,746,747,748,};
  static int qb_7_4[] = { 749,750,751,752,753,754,755,756,757,758,759,760,761,762,763,764,765,766,};
  static int qb_7_5[] = { 767,768,769,770,771,772,773,774,775,776,777,778,779,780,781,782,783,784,785,786,787,788,789,790,};
  static int qb_7_6[] = { 791,792,793,794,795,796,797,798,799,800,801,802,803,804,805,806,807,808,809,810,811,812,813,814,815,816,817,818,819,820,};
  static int qb_7_7[] = { 821,822,823,824,825,826,827,828,829,830,831,832,833,834,835,836,837,838,839,840,841,842,843,844,845,846,847,848,849,850,851,852,853,854,855,856,};
  static int qb_7_8[] = { 857,858,859,860,861,862,863,864,865,866,867,868,869,870,871,872,873,874,875,876,877,878,879,880,881,882,883,884,885,886,887,888,889,890,891,892,893,894,895,896,897,898,};
  static int qb_7_9[] = { 899,900,901,902,903,904,905,906,907,908,909,910,911,912,913,914,915,916,917,918,919,920,921,922,923,924,925,926,927,928,929,930,931,932,933,934,935,936,937,938,939,940,941,942,943,944,945,946,};
  static int qb_7_10[] = { 947,948,949,950,951,952,953,954,955,956,957,958,959,960,961,962,963,964,965,966,967,968,969,970,971,972,973,974,975,976,977,978,979,980,981,982,983,984,985,986,987,988,989,990,991,992,993,994,995,996,997,998,999,1000,};

  static int qb_8_2[] = { 1001,1002,1003,1004,1005,1006,1007,};
  static int qb_8_3[] = { 1008,1009,1010,1011,1012,1013,1014,1015,1016,1017,1018,1019,1020,1021,};
  static int qb_8_4[] = { 1022,1023,1024,1025,1026,1027,1028,1029,1030,1031,1032,1033,1034,1035,1036,1037,1038,1039,1040,1041,1042,};
  static int qb_8_5[] = { 1043,1044,1045,1046,1047,1048,1049,1050,1051,1052,1053,1054,1055,1056,1057,1058,1059,1060,1061,1062,1063,1064,1065,1066,1067,1068,1069,1070,};
  static int qb_8_6[] = { 1071,1072,1073,1074,1075,1076,1077,1078,1079,1080,1081,1082,1083,1084,1085,1086,1087,1088,1089,1090,1091,1092,1093,1094,1095,1096,1097,1098,1099,1100,1101,1102,1103,1104,1105,};
  static int qb_8_7[] = { 1106,1107,1108,1109,1110,1111,1112,1113,1114,1115,1116,1117,1118,1119,1120,1121,1122,1123,1124,1125,1126,1127,1128,1129,1130,1131,1132,1133,1134,1135,1136,1137,1138,1139,1140,1141,1142,1143,1144,1145,1146,1147,};
  static int qb_8_8[] = { 1148,1149,1150,1151,1152,1153,1154,1155,1156,1157,1158,1159,1160,1161,1162,1163,1164,1165,1166,1167,1168,1169,1170,1171,1172,1173,1174,1175,1176,1177,1178,1179,1180,1181,1182,1183,1184,1185,1186,1187,1188,1189,1190,1191,1192,1193,1194,1195,1196,};
  static int qb_8_9[] = { 1197,1198,1199,1200,1201,1202,1203,1204,1205,1206,1207,1208,1209,1210,1211,1212,1213,1214,1215,1216,1217,1218,1219,1220,1221,1222,1223,1224,1225,1226,1227,1228,1229,1230,1231,1232,1233,1234,1235,1236,1237,1238,1239,1240,1241,1242,1243,1244,1245,1246,1247,1248,1249,1250,1251,1252,};
  static int qb_8_10[] = { 1253,1254,1255,1256,1257,1258,1259,1260,1261,1262,1263,1264,1265,1266,1267,1268,1269,1270,1271,1272,1273,1274,1275,1276,1277,1278,1279,1280,1281,1282,1283,1284,1285,1286,1287,1288,1289,1290,1291,1292,1293,1294,1295,1296,1297,1298,1299,1300,1301,1302,1303,1304,1305,1306,1307,1308,1309,1310,1311,1312,1313,1314,1315,};

  static int qb_9_2[] = { 1316,1317,1318,1319,1320,1321,1322,1323,};
  static int qb_9_3[] = { 1324,1325,1326,1327,1328,1329,1330,1331,1332,1333,1334,1335,1336,1337,1338,1339,};
  static int qb_9_4[] = { 1340,1341,1342,1343,1344,1345,1346,1347,1348,1349,1350,1351,1352,1353,1354,1355,1356,1357,1358,1359,1360,1361,1362,1363,};
  static int qb_9_5[] = { 1364,1365,1366,1367,1368,1369,1370,1371,1372,1373,1374,1375,1376,1377,1378,1379,1380,1381,1382,1383,1384,1385,1386,1387,1388,1389,1390,1391,1392,1393,1394,1395,};
  static int qb_9_6[] = { 1396,1397,1398,1399,1400,1401,1402,1403,1404,1405,1406,1407,1408,1409,1410,1411,1412,1413,1414,1415,1416,1417,1418,1419,1420,1421,1422,1423,1424,1425,1426,1427,1428,1429,1430,1431,1432,1433,1434,1435,};
  static int qb_9_7[] = { 1436,1437,1438,1439,1440,1441,1442,1443,1444,1445,1446,1447,1448,1449,1450,1451,1452,1453,1454,1455,1456,1457,1458,1459,1460,1461,1462,1463,1464,1465,1466,1467,1468,1469,1470,1471,1472,1473,1474,1475,1476,1477,1478,1479,1480,1481,1482,1483,};
  static int qb_9_8[] = { 1484,1485,1486,1487,1488,1489,1490,1491,1492,1493,1494,1495,1496,1497,1498,1499,1500,1501,1502,1503,1504,1505,1506,1507,1508,1509,1510,1511,1512,1513,1514,1515,1516,1517,1518,1519,1520,1521,1522,1523,1524,1525,1526,1527,1528,1529,1530,1531,1532,1533,1534,1535,1536,1537,1538,1539,};
  static int qb_9_9[] = { 1540,1541,1542,1543,1544,1545,1546,1547,1548,1549,1550,1551,1552,1553,1554,1555,1556,1557,1558,1559,1560,1561,1562,1563,1564,1565,1566,1567,1568,1569,1570,1571,1572,1573,1574,1575,1576,1577,1578,1579,1580,1581,1582,1583,1584,1585,1586,1587,1588,1589,1590,1591,1592,1593,1594,1595,1596,1597,1598,1599,1600,1601,1602,1603,};
  static int qb_9_10[] = { 1604,1605,1606,1607,1608,1609,1610,1611,1612,1613,1614,1615,1616,1617,1618,1619,1620,1621,1622,1623,1624,1625,1626,1627,1628,1629,1630,1631,1632,1633,1634,1635,1636,1637,1638,1639,1640,1641,1642,1643,1644,1645,1646,1647,1648,1649,1650,1651,1652,1653,1654,1655,1656,1657,1658,1659,1660,1661,1662,1663,1664,1665,1666,1667,1668,1669,1670,1671,1672,1673,1674,1675,};

  static int qb_10_2[] = { 1676,1677,1678,1679,1680,1681,1682,1683,1684,};
  static int qb_10_3[] = { 1685,1686,1687,1688,1689,1690,1691,1692,1693,1694,1695,1696,1697,1698,1699,1700,1701,1702,};
  static int qb_10_4[] = { 1703,1704,1705,1706,1707,1708,1709,1710,1711,1712,1713,1714,1715,1716,1717,1718,1719,1720,1721,1722,1723,1724,1725,1726,1727,1728,1729,};
  static int qb_10_5[] = { 1730,1731,1732,1733,1734,1735,1736,1737,1738,1739,1740,1741,1742,1743,1744,1745,1746,1747,1748,1749,1750,1751,1752,1753,1754,1755,1756,1757,1758,1759,1760,1761,1762,1763,1764,1765,};
  static int qb_10_6[] = { 1766,1767,1768,1769,1770,1771,1772,1773,1774,1775,1776,1777,1778,1779,1780,1781,1782,1783,1784,1785,1786,1787,1788,1789,1790,1791,1792,1793,1794,1795,1796,1797,1798,1799,1800,1801,1802,1803,1804,1805,1806,1807,1808,1809,1810,};
  static int qb_10_7[] = { 1811,1812,1813,1814,1815,1816,1817,1818,1819,1820,1821,1822,1823,1824,1825,1826,1827,1828,1829,1830,1831,1832,1833,1834,1835,1836,1837,1838,1839,1840,1841,1842,1843,1844,1845,1846,1847,1848,1849,1850,1851,1852,1853,1854,1855,1856,1857,1858,1859,1860,1861,1862,1863,1864,};
  static int qb_10_8[] = { 1865,1866,1867,1868,1869,1870,1871,1872,1873,1874,1875,1876,1877,1878,1879,1880,1881,1882,1883,1884,1885,1886,1887,1888,1889,1890,1891,1892,1893,1894,1895,1896,1897,1898,1899,1900,1901,1902,1903,1904,1905,1906,1907,1908,1909,1910,1911,1912,1913,1914,1915,1916,1917,1918,1919,1920,1921,1922,1923,1924,1925,1926,1927,};
  static int qb_10_9[] = { 1928,1929,1930,1931,1932,1933,1934,1935,1936,1937,1938,1939,1940,1941,1942,1943,1944,1945,1946,1947,1948,1949,1950,1951,1952,1953,1954,1955,1956,1957,1958,1959,1960,1961,1962,1963,1964,1965,1966,1967,1968,1969,1970,1971,1972,1973,1974,1975,1976,1977,1978,1979,1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,};
  static int qb_10_10[] = { 2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,2023,2024,2025,2026,2027,2028,2029,2030,2031,2032,2033,2034,2035,2036,2037,2038,2039,2040,2041,2042,2043,2044,2045,2046,2047,2048,2049,2050,2051,2052,2053,2054,2055,2056,2057,2058,2059,2060,2061,2062,2063,2064,2065,2066,2067,2068,2069,2070,2071,2072,2073,2074,2075,2076,2077,2078,2079,2080,};




#define nullptr16 nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,  nullptr, nullptr, nullptr, nullptr, nullptr

int* eigen_quad_bubble_indices[] =
{
  nullptr, nullptr, nullptr,    nullptr,    nullptr,    nullptr,    nullptr,    nullptr,    nullptr,    nullptr,    nullptr,     nullptr, nullptr, nullptr, nullptr, nullptr, nullptr16,
  nullptr, nullptr, nullptr,    nullptr,    nullptr,    nullptr,    nullptr,    nullptr,    nullptr,    nullptr,    nullptr,     nullptr, nullptr, nullptr, nullptr, nullptr, nullptr16,
  nullptr, nullptr, qb_2_2,  qb_2_3,  qb_2_4,  qb_2_5,  qb_2_6,  qb_2_7,  qb_2_8,  qb_2_9,  qb_2_10,  nullptr, nullptr, nullptr, nullptr, nullptr, nullptr16,
  nullptr, nullptr, qb_3_2,  qb_3_3,  qb_3_4,  qb_3_5,  qb_3_6,  qb_3_7,  qb_3_8,  qb_3_9,  qb_3_10,  nullptr, nullptr, nullptr, nullptr, nullptr, nullptr16,
  nullptr, nullptr, qb_4_2,  qb_4_3,  qb_4_4,  qb_4_5,  qb_4_6,  qb_4_7,  qb_4_8,  qb_4_9,  qb_4_10,  nullptr, nullptr, nullptr, nullptr, nullptr, nullptr16,
  nullptr, nullptr, qb_5_2,  qb_5_3,  qb_5_4,  qb_5_5,  qb_5_6,  qb_5_7,  qb_5_8,  qb_5_9,  qb_5_10,  nullptr, nullptr, nullptr, nullptr, nullptr, nullptr16,
  nullptr, nullptr, qb_6_2,  qb_6_3,  qb_6_4,  qb_6_5,  qb_6_6,  qb_6_7,  qb_6_8,  qb_6_9,  qb_6_10,  nullptr, nullptr, nullptr, nullptr, nullptr, nullptr16,
  nullptr, nullptr, qb_7_2,  qb_7_3,  qb_7_4,  qb_7_5,  qb_7_6,  qb_7_7,  qb_7_8,  qb_7_9,  qb_7_10,  nullptr, nullptr, nullptr, nullptr, nullptr, nullptr16,
  nullptr, nullptr, qb_8_2,  qb_8_3,  qb_8_4,  qb_8_5,  qb_8_6,  qb_8_7,  qb_8_8,  qb_8_9,  qb_8_10,  nullptr, nullptr, nullptr, nullptr, nullptr, nullptr16,
  nullptr, nullptr, qb_9_2,  qb_9_3,  qb_9_4,  qb_9_5,  qb_9_6,  qb_9_7,  qb_9_8,  qb_9_9,  qb_9_10,  nullptr, nullptr, nullptr, nullptr, nullptr, nullptr16,
  nullptr, nullptr, qb_10_2, qb_10_3, qb_10_4, qb_10_5, qb_10_6, qb_10_7, qb_10_8, qb_10_9, qb_10_10, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr16
};


#define zero16  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

int eigen_quad_bubble_count[] =
{
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, zero16,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, zero16,
  0,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  0,  0,  0,  0,  0, zero16,
  0,  0,  2,  4,  6,  8, 10, 12, 14, 16, 18,  0,  0,  0,  0,  0, zero16,
  0,  0,  3,  6,  9, 12, 15, 18, 21, 24, 27,  0,  0,  0,  0,  0, zero16,
  0,  0,  4,  8, 12, 16, 20, 24, 28, 32, 36,  0,  0,  0,  0,  0, zero16,
  0,  0,  5, 10, 15, 20, 25, 30, 35, 40, 45,  0,  0,  0,  0,  0, zero16,
  0,  0,  6, 12, 18, 24, 30, 36, 42, 48, 54,  0,  0,  0,  0,  0, zero16,
  0,  0,  7, 14, 21, 28, 35, 42, 49, 56, 63,  0,  0,  0,  0,  0, zero16,
  0,  0,  8, 16, 24, 32, 40, 48, 56, 64, 72,  0,  0,  0,  0,  0, zero16,
  0,  0,  9, 18, 27, 36, 45, 54, 63, 72, 81,  0,  0,  0,  0,  0, zero16
};


int eigen_quad_vertex_indices[4] = { 0, 15, 16, 1 };

static int eigen_quad_edge_indices_0[22] =  {  0, 15, 15, 0, 30, 30, 32, 33, 36, 36, 38, 39, 42, 42, 44, 45, 48, 48, 50, 51, 54, 54 };
static int eigen_quad_edge_indices_1[22] =  { 15, 16, 16, 15, 17,17, 18, 19, 20, 20, 21, 22, 23, 23, 24, 25, 26, 26, 27, 28, 29, 29 };
static int eigen_quad_edge_indices_2[22] =  {  1, 16, 16, 1, 31, 31, 34, 35, 37, 37, 40, 41, 43, 43, 46, 47, 49, 49, 52, 53, 55, 55 };
static int eigen_quad_edge_indices_3[22] =  {  0,  1,  1, 0,  2,  2,  3,  4,  5,  5,  6,  7,  8,  8,  9, 10, 11, 11, 12, 13, 14, 14 };

int* eigen_quad_edge_indices[4] =
{
  eigen_quad_edge_indices_0,
  eigen_quad_edge_indices_1,
  eigen_quad_edge_indices_2,
  eigen_quad_edge_indices_3
};


#define oo H2D_MAKE_QUAD_ORDER
#define XX(a,b) oo(a,b), oo(a,b)

int eigen_quad_index_to_order[] =
{
  oo(1,1),  oo(1,1),  oo(1,2),  XX(1,3),  oo(1,4),  XX(1,5),  oo(1,6),  XX(1,7),  oo(1,8),  XX(1,9),  oo(1,10),
  oo(1,1),  oo(1,1),  oo(1,2),  XX(1,3),  oo(1,4),  XX(1,5),  oo(1,6),  XX(1,7),  oo(1,8),  XX(1,9),  oo(1,10),
  oo(2,1),  oo(2,1),  XX(3,1),  XX(3,1),  oo(4,1),  oo(4,1),  XX(5,1),  XX(5,1),  oo(6,1),  oo(6,1),  XX(7,1),
  XX(7,1),  oo(8,1),  oo(8,1),  XX(9,1),  XX(9,1),  oo(10,1), oo(10,1),

  oo(2,2), oo(2,3), oo(2,3), oo(2,4), oo(2,4), oo(2,4), oo(2,5), oo(2,5), oo(2,5), oo(2,5), oo(2,6), oo(2,6), oo(2,6), oo(2,6), oo(2,6), oo(2,7), oo(2,7), oo(2,7), oo(2,7), oo(2,7), oo(2,7), oo(2,8), oo(2,8), oo(2,8), oo(2,8), oo(2,8), oo(2,8), oo(2,8), oo(2,9), oo(2,9), oo(2,9), oo(2,9), oo(2,9), oo(2,9), oo(2,9), oo(2,9), oo(2,10), oo(2,10), oo(2,10), oo(2,10), oo(2,10), oo(2,10), oo(2,10), oo(2,10), oo(2,10),
  oo(3,2), oo(3,2), oo(3,3), oo(3,3), oo(3,3), oo(3,3), oo(3,4), oo(3,4), oo(3,4), oo(3,4), oo(3,4), oo(3,4), oo(3,5), oo(3,5), oo(3,5), oo(3,5), oo(3,5), oo(3,5), oo(3,5), oo(3,5), oo(3,6), oo(3,6), oo(3,6), oo(3,6), oo(3,6), oo(3,6), oo(3,6), oo(3,6), oo(3,6), oo(3,6), oo(3,7), oo(3,7), oo(3,7), oo(3,7), oo(3,7), oo(3,7), oo(3,7), oo(3,7), oo(3,7), oo(3,7), oo(3,7), oo(3,7), oo(3,8), oo(3,8), oo(3,8), oo(3,8), oo(3,8), oo(3,8), oo(3,8), oo(3,8), oo(3,8), oo(3,8), oo(3,8), oo(3,8), oo(3,8), oo(3,8), oo(3,9), oo(3,9), oo(3,9), oo(3,9), oo(3,9), oo(3,9), oo(3,9), oo(3,9), oo(3,9), oo(3,9), oo(3,9), oo(3,9), oo(3,9), oo(3,9), oo(3,9), oo(3,9), oo(3,10), oo(3,10), oo(3,10), oo(3,10), oo(3,10), oo(3,10), oo(3,10), oo(3,10), oo(3,10), oo(3,10), oo(3,10), oo(3,10), oo(3,10), oo(3,10), oo(3,10), oo(3,10), oo(3,10), oo(3,10),
  oo(4,2), oo(4,2), oo(4,2), oo(4,3), oo(4,3), oo(4,3), oo(4,3), oo(4,3), oo(4,3), oo(4,4), oo(4,4), oo(4,4), oo(4,4), oo(4,4), oo(4,4), oo(4,4), oo(4,4), oo(4,4), oo(4,5), oo(4,5), oo(4,5), oo(4,5), oo(4,5), oo(4,5), oo(4,5), oo(4,5), oo(4,5), oo(4,5), oo(4,5), oo(4,5), oo(4,6), oo(4,6), oo(4,6), oo(4,6), oo(4,6), oo(4,6), oo(4,6), oo(4,6), oo(4,6), oo(4,6), oo(4,6), oo(4,6), oo(4,6), oo(4,6), oo(4,6), oo(4,7), oo(4,7), oo(4,7), oo(4,7), oo(4,7), oo(4,7), oo(4,7), oo(4,7), oo(4,7), oo(4,7), oo(4,7), oo(4,7), oo(4,7), oo(4,7), oo(4,7), oo(4,7), oo(4,7), oo(4,7), oo(4,8), oo(4,8), oo(4,8), oo(4,8), oo(4,8), oo(4,8), oo(4,8), oo(4,8), oo(4,8), oo(4,8), oo(4,8), oo(4,8), oo(4,8), oo(4,8), oo(4,8), oo(4,8), oo(4,8), oo(4,8), oo(4,8), oo(4,8), oo(4,8), oo(4,9), oo(4,9), oo(4,9), oo(4,9), oo(4,9), oo(4,9), oo(4,9), oo(4,9), oo(4,9), oo(4,9), oo(4,9), oo(4,9), oo(4,9), oo(4,9), oo(4,9), oo(4,9), oo(4,9), oo(4,9), oo(4,9), oo(4,9), oo(4,9), oo(4,9), oo(4,9), oo(4,9), oo(4,10), oo(4,10), oo(4,10), oo(4,10), oo(4,10), oo(4,10), oo(4,10), oo(4,10), oo(4,10), oo(4,10), oo(4,10), oo(4,10), oo(4,10), oo(4,10), oo(4,10), oo(4,10), oo(4,10), oo(4,10), oo(4,10), oo(4,10), oo(4,10), oo(4,10), oo(4,10), oo(4,10), oo(4,10), oo(4,10), oo(4,10),
  oo(5,2), oo(5,2), oo(5,2), oo(5,2), oo(5,3), oo(5,3), oo(5,3), oo(5,3), oo(5,3), oo(5,3), oo(5,3), oo(5,3), oo(5,4), oo(5,4), oo(5,4), oo(5,4), oo(5,4), oo(5,4), oo(5,4), oo(5,4), oo(5,4), oo(5,4), oo(5,4), oo(5,4), oo(5,5), oo(5,5), oo(5,5), oo(5,5), oo(5,5), oo(5,5), oo(5,5), oo(5,5), oo(5,5), oo(5,5), oo(5,5), oo(5,5), oo(5,5), oo(5,5), oo(5,5), oo(5,5), oo(5,6), oo(5,6), oo(5,6), oo(5,6), oo(5,6), oo(5,6), oo(5,6), oo(5,6), oo(5,6), oo(5,6), oo(5,6), oo(5,6), oo(5,6), oo(5,6), oo(5,6), oo(5,6), oo(5,6), oo(5,6), oo(5,6), oo(5,6), oo(5,7), oo(5,7), oo(5,7), oo(5,7), oo(5,7), oo(5,7), oo(5,7), oo(5,7), oo(5,7), oo(5,7), oo(5,7), oo(5,7), oo(5,7), oo(5,7), oo(5,7), oo(5,7), oo(5,7), oo(5,7), oo(5,7), oo(5,7), oo(5,7), oo(5,7), oo(5,7), oo(5,7), oo(5,8), oo(5,8), oo(5,8), oo(5,8), oo(5,8), oo(5,8), oo(5,8), oo(5,8), oo(5,8), oo(5,8), oo(5,8), oo(5,8), oo(5,8), oo(5,8), oo(5,8), oo(5,8), oo(5,8), oo(5,8), oo(5,8), oo(5,8), oo(5,8), oo(5,8), oo(5,8), oo(5,8), oo(5,8), oo(5,8), oo(5,8), oo(5,8), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,9), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10), oo(5,10),
  oo(6,2), oo(6,2), oo(6,2), oo(6,2), oo(6,2), oo(6,3), oo(6,3), oo(6,3), oo(6,3), oo(6,3), oo(6,3), oo(6,3), oo(6,3), oo(6,3), oo(6,3), oo(6,4), oo(6,4), oo(6,4), oo(6,4), oo(6,4), oo(6,4), oo(6,4), oo(6,4), oo(6,4), oo(6,4), oo(6,4), oo(6,4), oo(6,4), oo(6,4), oo(6,4), oo(6,5), oo(6,5), oo(6,5), oo(6,5), oo(6,5), oo(6,5), oo(6,5), oo(6,5), oo(6,5), oo(6,5), oo(6,5), oo(6,5), oo(6,5), oo(6,5), oo(6,5), oo(6,5), oo(6,5), oo(6,5), oo(6,5), oo(6,5), oo(6,6), oo(6,6), oo(6,6), oo(6,6), oo(6,6), oo(6,6), oo(6,6), oo(6,6), oo(6,6), oo(6,6), oo(6,6), oo(6,6), oo(6,6), oo(6,6), oo(6,6), oo(6,6), oo(6,6), oo(6,6), oo(6,6), oo(6,6), oo(6,6), oo(6,6), oo(6,6), oo(6,6), oo(6,6), oo(6,7), oo(6,7), oo(6,7), oo(6,7), oo(6,7), oo(6,7), oo(6,7), oo(6,7), oo(6,7), oo(6,7), oo(6,7), oo(6,7), oo(6,7), oo(6,7), oo(6,7), oo(6,7), oo(6,7), oo(6,7), oo(6,7), oo(6,7), oo(6,7), oo(6,7), oo(6,7), oo(6,7), oo(6,7), oo(6,7), oo(6,7), oo(6,7), oo(6,7), oo(6,7), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,8), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,9), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10), oo(6,10),
  oo(7,2), oo(7,2), oo(7,2), oo(7,2), oo(7,2), oo(7,2), oo(7,3), oo(7,3), oo(7,3), oo(7,3), oo(7,3), oo(7,3), oo(7,3), oo(7,3), oo(7,3), oo(7,3), oo(7,3), oo(7,3), oo(7,4), oo(7,4), oo(7,4), oo(7,4), oo(7,4), oo(7,4), oo(7,4), oo(7,4), oo(7,4), oo(7,4), oo(7,4), oo(7,4), oo(7,4), oo(7,4), oo(7,4), oo(7,4), oo(7,4), oo(7,4), oo(7,5), oo(7,5), oo(7,5), oo(7,5), oo(7,5), oo(7,5), oo(7,5), oo(7,5), oo(7,5), oo(7,5), oo(7,5), oo(7,5), oo(7,5), oo(7,5), oo(7,5), oo(7,5), oo(7,5), oo(7,5), oo(7,5), oo(7,5), oo(7,5), oo(7,5), oo(7,5), oo(7,5), oo(7,6), oo(7,6), oo(7,6), oo(7,6), oo(7,6), oo(7,6), oo(7,6), oo(7,6), oo(7,6), oo(7,6), oo(7,6), oo(7,6), oo(7,6), oo(7,6), oo(7,6), oo(7,6), oo(7,6), oo(7,6), oo(7,6), oo(7,6), oo(7,6), oo(7,6), oo(7,6), oo(7,6), oo(7,6), oo(7,6), oo(7,6), oo(7,6), oo(7,6), oo(7,6), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,7), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,8), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,9), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10), oo(7,10),
  oo(8,2), oo(8,2), oo(8,2), oo(8,2), oo(8,2), oo(8,2), oo(8,2), oo(8,3), oo(8,3), oo(8,3), oo(8,3), oo(8,3), oo(8,3), oo(8,3), oo(8,3), oo(8,3), oo(8,3), oo(8,3), oo(8,3), oo(8,3), oo(8,3), oo(8,4), oo(8,4), oo(8,4), oo(8,4), oo(8,4), oo(8,4), oo(8,4), oo(8,4), oo(8,4), oo(8,4), oo(8,4), oo(8,4), oo(8,4), oo(8,4), oo(8,4), oo(8,4), oo(8,4), oo(8,4), oo(8,4), oo(8,4), oo(8,4), oo(8,5), oo(8,5), oo(8,5), oo(8,5), oo(8,5), oo(8,5), oo(8,5), oo(8,5), oo(8,5), oo(8,5), oo(8,5), oo(8,5), oo(8,5), oo(8,5), oo(8,5), oo(8,5), oo(8,5), oo(8,5), oo(8,5), oo(8,5), oo(8,5), oo(8,5), oo(8,5), oo(8,5), oo(8,5), oo(8,5), oo(8,5), oo(8,5), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,6), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,7), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,8), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,9), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10), oo(8,10),
  oo(9,2), oo(9,2), oo(9,2), oo(9,2), oo(9,2), oo(9,2), oo(9,2), oo(9,2), oo(9,3), oo(9,3), oo(9,3), oo(9,3), oo(9,3), oo(9,3), oo(9,3), oo(9,3), oo(9,3), oo(9,3), oo(9,3), oo(9,3), oo(9,3), oo(9,3), oo(9,3), oo(9,3), oo(9,4), oo(9,4), oo(9,4), oo(9,4), oo(9,4), oo(9,4), oo(9,4), oo(9,4), oo(9,4), oo(9,4), oo(9,4), oo(9,4), oo(9,4), oo(9,4), oo(9,4), oo(9,4), oo(9,4), oo(9,4), oo(9,4), oo(9,4), oo(9,4), oo(9,4), oo(9,4), oo(9,4), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,5), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,6), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,7), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,8), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,9), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10), oo(9,10),
  oo(10,2), oo(10,2), oo(10,2), oo(10,2), oo(10,2), oo(10,2), oo(10,2), oo(10,2), oo(10,2), oo(10,3), oo(10,3), oo(10,3), oo(10,3), oo(10,3), oo(10,3), oo(10,3), oo(10,3), oo(10,3), oo(10,3), oo(10,3), oo(10,3), oo(10,3), oo(10,3), oo(10,3), oo(10,3), oo(10,3), oo(10,3), oo(10,4), oo(10,4), oo(10,4), oo(10,4), oo(10,4), oo(10,4), oo(10,4), oo(10,4), oo(10,4), oo(10,4), oo(10,4), oo(10,4), oo(10,4), oo(10,4), oo(10,4), oo(10,4), oo(10,4), oo(10,4), oo(10,4), oo(10,4), oo(10,4), oo(10,4), oo(10,4), oo(10,4), oo(10,4), oo(10,4), oo(10,4), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,5), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,6), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,7), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,8), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,9), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10), oo(10,10),
};


static Shapeset::shape_fn_t* eigen_quad_shape_fn_table[1] =
{
  eigen_quad_fn
};

static Shapeset::shape_fn_t* eigen_quad_shape_fn_table_x[1] =
{
  eigen_quad_fn_dx
};

static Shapeset::shape_fn_t* eigen_quad_shape_fn_table_y[1] =
{
  eigen_quad_fn_dy
};

//// triangle tables and class constructor ///////////////////////////////////////////////

static Shapeset::shape_fn_t** eigen_shape_fn_table[2] =
{
  eigen_quad_shape_fn_table,
  eigen_quad_shape_fn_table
};


static Shapeset::shape_fn_t** eigen_shape_fn_table_x[2] =
{
  nullptr,
  eigen_quad_shape_fn_table_x
};

static Shapeset::shape_fn_t** eigen_shape_fn_table_y[2] =
{
  nullptr,
  eigen_quad_shape_fn_table_y
};

static int* eigen_vertex_indices[2] =
{
  nullptr,
  eigen_quad_vertex_indices
};


static int** eigen_edge_indices[2] =
{
  eigen_quad_edge_indices,
  eigen_quad_edge_indices
};

/*static int* eigen_edge_count[2] =
{
  nullptr,
  eigen_quad_edge_count
};
*/

static int** eigen_bubble_indices[2] =
{
  nullptr,
  eigen_quad_bubble_indices
};


static int* eigen_bubble_count[2] =
{
  nullptr,
  eigen_quad_bubble_count
};


static int* eigen_index_to_order[2] =
{
  nullptr,
  eigen_quad_index_to_order
};


H1ShapesetEigen::H1ShapesetEigen()
{
  shape_table[0] = eigen_shape_fn_table;
  shape_table[1] = eigen_shape_fn_table_x;
  shape_table[2] = eigen_shape_fn_table_y;
  shape_table[3] = nullptr;
  shape_table[4] = nullptr;
  shape_table[5] = nullptr;

  vertex_indices = eigen_vertex_indices;
  edge_indices = eigen_edge_indices;
  bubble_indices = eigen_bubble_indices;
  bubble_count = eigen_bubble_count;
  index_to_order = eigen_index_to_order;

  ref_vert[0][0][0] = -1.0;
  ref_vert[0][0][1] = -1.0;
  ref_vert[0][1][0] =  1.0;
  ref_vert[0][1][1] = -1.0;
  ref_vert[0][2][0] = -1.0;
  ref_vert[0][2][1] =  1.0;

  ref_vert[1][0][0] = -1.0;
  ref_vert[1][0][1] = -1.0;
  ref_vert[1][1][0] =  1.0;
  ref_vert[1][1][1] = -1.0;
  ref_vert[1][2][0] =  1.0;
  ref_vert[1][2][1] =  1.0;
  ref_vert[1][3][0] = -1.0;
  ref_vert[1][3][1] =  1.0;

  max_order = 10;
  num_components = 1;

  max_index[0] = 2080;
  max_index[1] = 2080;

  cei[0] = 27;
  cei[1] = 3;
  cei[2] = -2;

  comb_table = nullptr;

  set_mode(HERMES_MODE_TRIANGLE);
}
