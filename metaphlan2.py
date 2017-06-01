#!/usr/bin/env python
from __future__ import with_statement 

# ==============================================================================
# MetaPhlAn v2.x: METAgenomic PHyLogenetic ANalysis for taxonomic classification
#                 of metagenomic data
#
# Authors: Nicola Segata (nicola.segata@unitn.it), 
#          Duy Tin Truong (duytin.truong@unitn.it)
#
# Please type "./metaphlan2.py -h" for usage help
#
# ==============================================================================

__author__ = 'Nicola Segata (nicola.segata@unitn.it), Duy Tin Truong (duytin.truong@unitn.it)'
__version__ = '2.6.0'
__date__ = '19 August 2016'


import sys
import os
import stat
import re
from binascii import b2a_uu 

try:
    import numpy as np 
except ImportError:
    sys.stderr.write("Error! numpy python library not detected!!\n")
    sys.exit(1)
import tempfile as tf
import argparse as ap
import subprocess as subp
import multiprocessing as mp
from collections import defaultdict as defdict
import bz2 
import itertools
from distutils.version import LooseVersion
try:
    import cPickle as pickle
except:
    import pickle

   
#**********************************************************************************************
#  Modification of Code :                                                                     *
#  Modified the code so instead of using the current clade IDs, which are numbers, we will    *
#      use the clade_names                                                                    *
#      Users reported the biom output is invalid and also the IDs were changing from run to   *
#      run.                                                                                   *
#  George Weingart    05/22/2017   george.weingart@mail.com                                   *
#**********************************************************************************************



#*************************************************************
#*  Imports related to biom file generation                  *
#*************************************************************
try:
    import biom
    import biom.table
    import numpy as np
except ImportError:
    sys.stderr.write("Warning! Biom python library not detected!"
                     "\n Exporting to biom format will not work!\n")
try:
    import json
except ImportError:
    sys.stderr.write("Warning! json python library not detected!"
                     "\n Exporting to biom format will not work!\n")

# This set contains the markers that after careful validation are found to have low precision or recall
# We esclude the markers here to avoid generating a new marker DB when changing just few markers
markers_to_exclude = \
    set([
        'NC_001782.1','GeneID:17099689','gi|419819595|ref|NZ_AJRE01000517.1|:1-118',
        'GeneID:10498696', 'GeneID:10498710', 'GeneID:10498726', 'GeneID:10498735',
        'GeneID:10498757', 'GeneID:10498760', 'GeneID:10498761', 'GeneID:10498763',
        'GeneID:11294465', 'GeneID:14181982', 'GeneID:14182132', 'GeneID:14182146',
        'GeneID:14182148', 'GeneID:14182328', 'GeneID:14182639', 'GeneID:14182647',
        'GeneID:14182650', 'GeneID:14182663', 'GeneID:14182683', 'GeneID:14182684',
        'GeneID:14182691', 'GeneID:14182803', 'GeneID:14296322', 'GeneID:1489077',
        'GeneID:1489080', 'GeneID:1489081', 'GeneID:1489084', 'GeneID:1489085',
        'GeneID:1489088', 'GeneID:1489089', 'GeneID:1489090', 'GeneID:1489528',
        'GeneID:1489530', 'GeneID:1489531', 'GeneID:1489735', 'GeneID:1491873',
        'GeneID:1491889', 'GeneID:1491962', 'GeneID:1491963', 'GeneID:1491964',
        'GeneID:1491965', 'GeneID:17099689', 'GeneID:1724732', 'GeneID:17494231',
        'GeneID:2546403', 'GeneID:2703374', 'GeneID:2703375', 'GeneID:2703498',
        'GeneID:2703531', 'GeneID:2772983', 'GeneID:2772989', 'GeneID:2772991',
        'GeneID:2772993', 'GeneID:2772995', 'GeneID:2773037', 'GeneID:2777387',
        'GeneID:2777399', 'GeneID:2777400', 'GeneID:2777439', 'GeneID:2777493',
        'GeneID:2777494', 'GeneID:3077424', 'GeneID:3160801', 'GeneID:3197323',
        'GeneID:3197355', 'GeneID:3197400', 'GeneID:3197428', 'GeneID:3783722',
        'GeneID:3783750', 'GeneID:3953004', 'GeneID:3959334', 'GeneID:3964368',
        'GeneID:3964370', 'GeneID:4961452', 'GeneID:5075645', 'GeneID:5075646',
        'GeneID:5075647', 'GeneID:5075648', 'GeneID:5075649', 'GeneID:5075650',
        'GeneID:5075651', 'GeneID:5075652', 'GeneID:5075653', 'GeneID:5075654',
        'GeneID:5075655', 'GeneID:5075656', 'GeneID:5075657', 'GeneID:5075658',
        'GeneID:5075659', 'GeneID:5075660', 'GeneID:5075661', 'GeneID:5075662',
        'GeneID:5075663', 'GeneID:5075664', 'GeneID:5075665', 'GeneID:5075667',
        'GeneID:5075668', 'GeneID:5075669', 'GeneID:5075670', 'GeneID:5075671',
        'GeneID:5075672', 'GeneID:5075673', 'GeneID:5075674', 'GeneID:5075675',
        'GeneID:5075676', 'GeneID:5075677', 'GeneID:5075678', 'GeneID:5075679',
        'GeneID:5075680', 'GeneID:5075681', 'GeneID:5075682', 'GeneID:5075683',
        'GeneID:5075684', 'GeneID:5075685', 'GeneID:5075686', 'GeneID:5075687',
        'GeneID:5075688', 'GeneID:5075689', 'GeneID:5075690', 'GeneID:5075691',
        'GeneID:5075692', 'GeneID:5075693', 'GeneID:5075694', 'GeneID:5075695',
        'GeneID:5075696', 'GeneID:5075697', 'GeneID:5075698', 'GeneID:5075700',
        'GeneID:5075701', 'GeneID:5075702', 'GeneID:5075703', 'GeneID:5075704',
        'GeneID:5075705', 'GeneID:5075707', 'GeneID:5075708', 'GeneID:5075709',
        'GeneID:5075710', 'GeneID:5075711', 'GeneID:5075712', 'GeneID:5075713',
        'GeneID:5075714', 'GeneID:5075715', 'GeneID:5075716', 'GeneID:5176189',
        'GeneID:6803896', 'GeneID:6803915', 'GeneID:7944151', 'GeneID:927334',
        'GeneID:927335', 'GeneID:927337', 'GeneID:940263', 'GeneID:9538324',
        'NC_003977.1', 'gi|103485498|ref|NC_008048.1|:1941166-1942314',
        'gi|108802856|ref|NC_008148.1|:1230231-1230875',
        'gi|124806686|ref|XM_001350760.1|',
        'gi|126661648|ref|NZ_AAXW01000149.1|:c1513-1341',
        'gi|149172845|ref|NZ_ABBW01000029.1|:970-1270',
        'gi|153883242|ref|NZ_ABDQ01000074.1|:79-541',
        'gi|167031021|ref|NC_010322.1|:1834668-1835168',
        'gi|171344510|ref|NZ_ABJO01001391.1|:1-116',
        'gi|171346813|ref|NZ_ABJO01001728.1|:c109-1',
        'gi|190640924|ref|NZ_ABRC01000948.1|:c226-44',
        'gi|223045343|ref|NZ_ACEN01000042.1|:1-336',
        'gi|224580998|ref|NZ_GG657387.1|:c114607-114002',
        'gi|224993759|ref|NZ_ACFY01000068.1|:c357-1',
        'gi|237784637|ref|NC_012704.1|:141000-142970',
        'gi|237784637|ref|NC_012704.1|:c2048315-2047083',
        'gi|240136783|ref|NC_012808.1|:1928224-1928961',
        'gi|255319020|ref|NZ_ACVR01000025.1|:28698-29132',
        'gi|260590341|ref|NZ_ACEO02000062.1|:c387-151',
        'gi|262368201|ref|NZ_GG704964.1|:733100-733978',
        'gi|262369811|ref|NZ_GG704966.1|:c264858-264520',
        'gi|288559258|ref|NC_013790.1|:448046-451354',
        'gi|288559258|ref|NC_013790.1|:532047-533942',
        'gi|294794157|ref|NZ_GG770200.1|:245344-245619',
        'gi|304372805|ref|NC_014448.1|:444677-445120',
        'gi|304372805|ref|NC_014448.1|:707516-708268',
        'gi|304372805|ref|NC_014448.1|:790263-792257',
        'gi|304372805|ref|NC_014448.1|:c367313-364470',
        'gi|304372805|ref|NC_014448.1|:c659144-658272',
        'gi|304372805|ref|NC_014448.1|:c772578-770410',
        'gi|304372805|ref|NC_014448.1|:c777901-777470',
        'gi|306477407|ref|NZ_GG770409.1|:c1643877-1643338',
        'gi|317120849|ref|NC_014831.1|:c891121-890144',
        'gi|323356441|ref|NZ_GL698442.1|:560-682',
        'gi|324996766|ref|NZ_BABV01000451.1|:10656-11579',
        'gi|326579405|ref|NZ_AEGQ01000006.1|:2997-3791',
        'gi|326579407|ref|NZ_AEGQ01000008.1|:c45210-44497',
        'gi|326579433|ref|NZ_AEGQ01000034.1|:346-3699',
        'gi|329889017|ref|NZ_GL883086.1|:586124-586804',
        'gi|330822653|ref|NC_015422.1|:2024431-2025018',
        'gi|335053104|ref|NZ_AFIL01000010.1|:c33862-32210',
        'gi|339304121|ref|NZ_AEOR01000258.1|:c294-1',
        'gi|339304277|ref|NZ_AEOR01000414.1|:1-812',
        'gi|342211239|ref|NZ_AFUK01000001.1|:790086-790835',
        'gi|342211239|ref|NZ_AFUK01000001.1|:c1579497-1578787',
        'gi|342213707|ref|NZ_AFUJ01000005.1|:48315-48908',
        'gi|355707189|ref|NZ_JH376566.1|:326756-326986',
        'gi|355707384|ref|NZ_JH376567.1|:90374-91453',
        'gi|355707384|ref|NZ_JH376567.1|:c388018-387605',
        'gi|355708440|ref|NZ_JH376569.1|:c80380-79448',
        'gi|358051729|ref|NZ_AEUN01000100.1|:c120-1',
        'gi|365983217|ref|XM_003668394.1|',
        'gi|377571722|ref|NZ_BAFD01000110.1|:c1267-29',
        'gi|377684864|ref|NZ_CM001194.1|:c1159954-1159619',
        'gi|377684864|ref|NZ_CM001194.1|:c4966-4196',
        'gi|378759497|ref|NZ_AFXE01000152.1|:1628-2215',
        'gi|378835506|ref|NC_016829.1|:112560-113342',
        'gi|378835506|ref|NC_016829.1|:114945-115193',
        'gi|378835506|ref|NC_016829.1|:126414-127151',
        'gi|378835506|ref|NC_016829.1|:272056-272403',
        'gi|378835506|ref|NC_016829.1|:272493-272786',
        'gi|378835506|ref|NC_016829.1|:358647-360863',
        'gi|378835506|ref|NC_016829.1|:37637-38185',
        'gi|378835506|ref|NC_016829.1|:60012-60497',
        'gi|378835506|ref|NC_016829.1|:606819-607427',
        'gi|378835506|ref|NC_016829.1|:607458-607760',
        'gi|378835506|ref|NC_016829.1|:826192-826821',
        'gi|378835506|ref|NC_016829.1|:c451932-451336',
        'gi|378835506|ref|NC_016829.1|:c460520-459951',
        'gi|378835506|ref|NC_016829.1|:c483843-482842',
        'gi|378835506|ref|NC_016829.1|:c544660-543638',
        'gi|378835506|ref|NC_016829.1|:c556383-555496',
        'gi|378835506|ref|NC_016829.1|:c632166-631228',
        'gi|378835506|ref|NC_016829.1|:c805066-802691',
        'gi|384124469|ref|NC_017160.1|:c2157447-2156863',
        'gi|385263288|ref|NZ_AJST01000001.1|:594143-594940',
        'gi|385858114|ref|NC_017519.1|:10252-10746',
        'gi|385858114|ref|NC_017519.1|:104630-104902',
        'gi|385858114|ref|NC_017519.1|:154292-156016',
        'gi|385858114|ref|NC_017519.1|:205158-206462',
        'gi|385858114|ref|NC_017519.1|:507239-507703',
        'gi|385858114|ref|NC_017519.1|:518924-519772',
        'gi|385858114|ref|NC_017519.1|:524712-525545',
        'gi|385858114|ref|NC_017519.1|:528387-528785',
        'gi|385858114|ref|NC_017519.1|:532275-533429',
        'gi|385858114|ref|NC_017519.1|:586402-586824',
        'gi|385858114|ref|NC_017519.1|:621696-622226',
        'gi|385858114|ref|NC_017519.1|:673673-676105',
        'gi|385858114|ref|NC_017519.1|:706602-708218',
        'gi|385858114|ref|NC_017519.1|:710627-711997',
        'gi|385858114|ref|NC_017519.1|:744974-745456',
        'gi|385858114|ref|NC_017519.1|:791055-791801',
        'gi|385858114|ref|NC_017519.1|:805643-807430',
        'gi|385858114|ref|NC_017519.1|:c172050-170809',
        'gi|385858114|ref|NC_017519.1|:c334545-333268',
        'gi|385858114|ref|NC_017519.1|:c383474-383202',
        'gi|385858114|ref|NC_017519.1|:c450880-450389',
        'gi|385858114|ref|NC_017519.1|:c451975-451001',
        'gi|385858114|ref|NC_017519.1|:c470488-470036',
        'gi|385858114|ref|NC_017519.1|:c485596-484598',
        'gi|385858114|ref|NC_017519.1|:c58658-58065',
        'gi|385858114|ref|NC_017519.1|:c592754-591081',
        'gi|385858114|ref|NC_017519.1|:c59590-58820',
        'gi|385858114|ref|NC_017519.1|:c601339-600575',
        'gi|385858114|ref|NC_017519.1|:c76080-75160',
        'gi|385858114|ref|NC_017519.1|:c97777-96302',
        'gi|391227518|ref|NZ_CM001514.1|:c1442504-1440237',
        'gi|391227518|ref|NZ_CM001514.1|:c3053472-3053023',
        'gi|394749766|ref|NZ_AHHC01000069.1|:3978-6176',
        'gi|398899615|ref|NZ_AKJK01000021.1|:28532-29209',
        'gi|406580057|ref|NZ_AJRD01000017.1|:c17130-15766',
        'gi|406584668|ref|NZ_AJQZ01000017.1|:c1397-771',
        'gi|408543458|ref|NZ_AJLO01000024.1|:67702-68304',
        'gi|410936685|ref|NZ_AJRF02000012.1|:21785-22696',
        'gi|41406098|ref|NC_002944.2|:c4468304-4467864',
        'gi|416998679|ref|NZ_AEXI01000003.1|:c562937-562176',
        'gi|417017738|ref|NZ_AEYL01000489.1|:c111-1',
        'gi|417018375|ref|NZ_AEYL01000508.1|:100-238',
        'gi|418576506|ref|NZ_AHKB01000025.1|:c7989-7669',
        'gi|419819595|ref|NZ_AJRE01000517.1|:1-118',
        'gi|421806549|ref|NZ_AMTB01000006.1|:c181247-180489',
        'gi|422320815|ref|NZ_GL636045.1|:28704-29048',
        'gi|422320874|ref|NZ_GL636046.1|:4984-5742',
        'gi|422323244|ref|NZ_GL636061.1|:479975-480520',
        'gi|422443048|ref|NZ_GL383112.1|:663738-664823',
        'gi|422552858|ref|NZ_GL383469.1|:c216727-215501',
        'gi|422859491|ref|NZ_GL878548.1|:c271832-271695',
        'gi|423012810|ref|NZ_GL982453.1|:3888672-3888935',
        'gi|423012810|ref|NZ_GL982453.1|:4541873-4542328',
        'gi|423012810|ref|NZ_GL982453.1|:c2189976-2188582',
        'gi|423012810|ref|NZ_GL982453.1|:c5471232-5470300',
        'gi|423262555|ref|NC_019552.1|:24703-25212',
        'gi|423262555|ref|NC_019552.1|:28306-30696',
        'gi|423262555|ref|NC_019552.1|:284252-284581',
        'gi|423262555|ref|NC_019552.1|:311161-311373',
        'gi|423262555|ref|NC_019552.1|:32707-34497',
        'gi|423262555|ref|NC_019552.1|:34497-35237',
        'gi|423262555|ref|NC_019552.1|:53691-56813',
        'gi|423262555|ref|NC_019552.1|:c388986-386611',
        'gi|423262555|ref|NC_019552.1|:c523106-522528',
        'gi|423689090|ref|NZ_CM001513.1|:c1700632-1699448',
        'gi|423689090|ref|NZ_CM001513.1|:c1701670-1700651',
        'gi|423689090|ref|NZ_CM001513.1|:c5739118-5738390',
        'gi|427395956|ref|NZ_JH992914.1|:c592682-591900',
        'gi|427407324|ref|NZ_JH992904.1|:c2681223-2679463',
        'gi|451952303|ref|NZ_AJRB03000021.1|:1041-1574',
        'gi|452231579|ref|NZ_AEKA01000123.1|:c18076-16676',
        'gi|459791914|ref|NZ_CM001824.1|:c899379-899239',
        'gi|471265562|ref|NC_020815.1|:3155799-3156695',
        'gi|472279780|ref|NZ_ALPV02000001.1|:33911-36751',
        'gi|482733945|ref|NZ_AHGZ01000071.1|:10408-11154',
        'gi|483051300|ref|NZ_ALYK02000034.1|:c37582-36650',
        'gi|483051300|ref|NZ_ALYK02000034.1|:c38037-37582',
        'gi|483993347|ref|NZ_AMXG01000045.1|:251724-253082',
        'gi|484100856|ref|NZ_JH670250.1|:600643-602949',
        'gi|484115941|ref|NZ_AJXG01000093.1|:567-947',
        'gi|484228609|ref|NZ_JH730929.1|:c103784-99021',
        'gi|484228797|ref|NZ_JH730960.1|:c16193-12429',
        'gi|484228814|ref|NZ_JH730962.1|:c29706-29260',
        'gi|484228929|ref|NZ_JH730981.1|:18645-22060',
        'gi|484228939|ref|NZ_JH730983.1|:42943-43860',
        'gi|484266598|ref|NZ_AKGC01000024.1|:118869-119636',
        'gi|484327375|ref|NZ_AKVP01000093.1|:1-1281',
        'gi|484328234|ref|NZ_AKVP01000127.1|:c325-110',
        'gi|487376144|ref|NZ_KB911257.1|:600445-601482',
        'gi|487376194|ref|NZ_KB911260.1|:146228-146533',
        'gi|487381776|ref|NZ_KB911485.1|:101242-103083',
        'gi|487381776|ref|NZ_KB911485.1|:c32472-31627',
        'gi|487381800|ref|NZ_KB911486.1|:39414-39872',
        'gi|487381828|ref|NZ_KB911487.1|:15689-17026',
        'gi|487381846|ref|NZ_KB911488.1|:13678-13821',
        'gi|487382089|ref|NZ_KB911497.1|:23810-26641',
        'gi|487382176|ref|NZ_KB911501.1|:c497-381',
        'gi|487382213|ref|NZ_KB911502.1|:12706-13119',
        'gi|487382247|ref|NZ_KB911505.1|:c7595-6663',
        'gi|490551798|ref|NZ_AORG01000011.1|:40110-41390',
        'gi|491099398|ref|NZ_KB849654.1|:c720460-719912',
        'gi|491124812|ref|NZ_KB849705.1|:1946500-1946937',
        'gi|491155563|ref|NZ_KB849732.1|:46469-46843',
        'gi|491155563|ref|NZ_KB849732.1|:46840-47181',
        'gi|491155563|ref|NZ_KB849732.1|:47165-48616',
        'gi|491155563|ref|NZ_KB849732.1|:55055-56662',
        'gi|491155563|ref|NZ_KB849732.1|:56662-57351',
        'gi|491155563|ref|NZ_KB849732.1|:6101-7588',
        'gi|491155563|ref|NZ_KB849732.1|:7657-8073',
        'gi|491349766|ref|NZ_KB850082.1|:441-941',
        'gi|491395079|ref|NZ_KB850142.1|:1461751-1462554',
        'gi|512608407|ref|NZ_KE150401.1|:c156891-156016',
        'gi|518653462|ref|NZ_ATLM01000004.1|:c89669-89247',
        'gi|520818261|ref|NZ_ATLQ01000015.1|:480744-481463',
        'gi|520822538|ref|NZ_ATLQ01000063.1|:103173-103283',
        'gi|520826510|ref|NZ_ATLQ01000092.1|:c13892-13563',
        'gi|544644736|ref|NZ_KE747865.1|:68388-69722',
        'gi|545347918|ref|NZ_KE952096.1|:c83873-81831',
        'gi|550735774|gb|AXMM01000002.1|:c743886-743575',
        'gi|552875787|ref|NZ_KI515684.1|:c584270-583890',
        'gi|552876418|ref|NZ_KI515685.1|:36713-37258',
        'gi|552876418|ref|NZ_KI515685.1|:432422-433465',
        'gi|552876418|ref|NZ_KI515685.1|:c1014617-1014117',
        'gi|552876418|ref|NZ_KI515685.1|:c931935-931327',
        'gi|552876815|ref|NZ_KI515686.1|:613740-614315',
        'gi|552879811|ref|NZ_AXME01000001.1|:1146402-1146932',
        'gi|552879811|ref|NZ_AXME01000001.1|:40840-41742',
        'gi|552879811|ref|NZ_AXME01000001.1|:49241-49654',
        'gi|552891898|ref|NZ_AXMG01000001.1|:99114-99290',
        'gi|552891898|ref|NZ_AXMG01000001.1|:c1460921-1460529',
        'gi|552895565|ref|NZ_AXMI01000001.1|:619555-620031',
        'gi|552895565|ref|NZ_AXMI01000001.1|:c14352-13837',
        'gi|552896371|ref|NZ_AXMI01000002.1|:c148595-146280',
        'gi|552897201|ref|NZ_AXMI01000004.1|:c231437-230883',
        'gi|552902020|ref|NZ_AXMK01000001.1|:c1625038-1624022',
        'gi|556346902|ref|NZ_KI535485.1|:c828278-827901',
        'gi|556478613|ref|NZ_KI535633.1|:3529392-3530162',
        'gi|560534311|ref|NZ_AYSF01000111.1|:26758-29049',
        'gi|564165687|gb|AYLX01000355.1|:10906-11166',
        'gi|564169776|gb|AYLX01000156.1|:1-185',
        'gi|564938696|gb|AWYH01000018.1|:c75674-75039', 'gi|67993724|ref|XM_664440.1|',
        'gi|68059117|ref|XM_666447.1|', 'gi|68062389|ref|XM_668109.1|',
        'gi|71730848|gb|AAAM03000019.1|:c14289-12877', 'gi|82753723|ref|XM_722699.1|',
        'gi|82775382|ref|NC_007606.1|:2249487-2250014', 'gi|82793634|ref|XM_723027.1|'
        ])

tax_units = "kpcofgst"

if float(sys.version_info[0]) < 3.0:
    def read_and_split( ofn  ):
        return (l.strip().split('\t') for l in ofn)
    def read_and_split_line( line ):
        return line.strip().split('\t')
else:
    def read_and_split( ofn ):
        return (str(l,encoding='utf-8').strip().split('\t') for l in ofn)
    def read_and_split_line( line ):
        return str(line,encoding='utf-8').strip().split('\t')


def plain_read_and_split( ofn ):
    return (l.strip().split('\t') for l in ofn)

def plain_read_and_split_line( l ):
    return l.strip().split('\t')



if float(sys.version_info[0]) < 3.0:
    def mybytes( val ):
        return val
else:
    def mybytes( val ):
        return bytes(val,encoding='utf-8')
    
# get the directory that contains this script
metaphlan2_script_install_folder=os.path.dirname(os.path.abspath(__file__))

def read_params(args):
    p = ap.ArgumentParser( description= 
            "DESCRIPTION\n"
            " MetaPhlAn version "+__version__+" ("+__date__+"): \n"
            " METAgenomic PHyLogenetic ANalysis for metagenomic taxonomic profiling.\n\n"
            "AUTHORS: "+__author__+"\n\n"
            "COMMON COMMANDS\n\n"
            " We assume here that metaphlan2.py is in the system path and that mpa_dir bash variable contains the\n"
            " main MetaPhlAn folder. Also BowTie2 should be in the system path with execution and read\n"
            " permissions, and Perl should be installed)\n\n"
           
            "\n========== MetaPhlAn 2 clade-abundance estimation ================= \n\n"
            "The basic usage of MetaPhlAn 2 consists in the identification of the clades (from phyla to species and \n"
            "strains in particular cases) present in the metagenome obtained from a microbiome sample and their \n"
            "relative abundance. This correspond to the default analysis type (--analysis_type rel_ab).\n\n"

            "*  Profiling a metagenome from raw reads:\n"
            "$ metaphlan2.py metagenome.fastq --input_type fastq\n\n"
            
            "*  You can take advantage of multiple CPUs and save the intermediate BowTie2 output for re-running\n"
            "   MetaPhlAn extremely quickly:\n"
            "$ metaphlan2.py metagenome.fastq --bowtie2out metagenome.bowtie2.bz2 --nproc 5 --input_type fastq\n\n"
            
            "*  If you already mapped your metagenome against the marker DB (using a previous MetaPhlAn run), you\n"
            "   can obtain the results in few seconds by using the previously saved --bowtie2out file and \n"
            "   specifying the input (--input_type bowtie2out):\n"
            "$ metaphlan2.py metagenome.bowtie2.bz2 --nproc 5 --input_type bowtie2out\n\n"
            
            "*  You can also provide an externally BowTie2-mapped SAM if you specify this format with \n"
            "   --input_type. Two steps: first apply BowTie2 and then feed MetaPhlAn2 with the obtained sam:\n"
            "$ bowtie2 --sam-no-hd --sam-no-sq --no-unal --very-sensitive -S metagenome.sam -x ${mpa_dir}/db_v20/mpa_v20_m200 -U metagenome.fastq\n"
            "$ metaphlan2.py metagenome.sam --input_type sam > profiled_metagenome.txt\n\n"
            
            "*  Multiple alternative ways to pass the input are also available:\n"
            "$ cat metagenome.fastq | metaphlan2.py --input_type fastq \n"
            "$ tar xjf metagenome.tar.bz2 --to-stdout | metaphlan2.py --input_type fastq \n"
            "$ metaphlan2.py --input_type fastq < metagenome.fastq\n"
            "$ metaphlan2.py --input_type fastq <(bzcat metagenome.fastq.bz2)\n"
            "$ metaphlan2.py --input_type fastq <(zcat metagenome_1.fastq.gz metagenome_2.fastq.gz)\n\n"

            "*  We can also natively handle paired-end metagenomes, and, more generally, metagenomes stored in \n"
            "  multiple files (but you need to specify the --bowtie2out parameter):\n"
            "$ metaphlan2.py metagenome_1.fastq,metagenome_2.fastq --bowtie2out metagenome.bowtie2.bz2 --nproc 5 --input_type fastq\n\n"
            "\n------------------------------------------------------------------- \n \n\n"
        
            
            "\n========== MetaPhlAn 2 strain tracking ============================ \n\n"
            "MetaPhlAn 2 introduces the capability of charachterizing organisms at the strain level using non\n"
            "aggregated marker information. Such capability comes with several slightly different flavours and \n"
            "are a way to perform strain tracking and comparison across multiple samples.\n"
            "Usually, MetaPhlAn 2 is first ran with the default --analysis_type to profile the species present in\n"
            "the community, and then a strain-level profiling can be performed to zoom-in into specific species\n"
            "of interest. This operation can be performed quickly as it exploits the --bowtie2out intermediate \n"
            "file saved during the execution of the default analysis type.\n\n"
           
            "*  The following command will output the abundance of each marker with a RPK (reads per kil-base) \n"
            "   higher 0.0. (we are assuming that metagenome_outfmt.bz2 has been generated before as \n"
            "   shown above).\n"
            "$ metaphlan2.py -t marker_ab_table metagenome_outfmt.bz2 --input_type bowtie2out > marker_abundance_table.txt\n"
            "   The obtained RPK can be optionally normalized by the total number of reads in the metagenome \n"
            "   to guarantee fair comparisons of abundances across samples. The number of reads in the metagenome\n"
            "   needs to be passed with the '--nreads' argument\n\n"

            "*  The list of markers present in the sample can be obtained with '-t marker_pres_table'\n"
            "$ metaphlan2.py -t marker_pres_table metagenome_outfmt.bz2 --input_type bowtie2out > marker_abundance_table.txt\n"
            "   The --pres_th argument (default 1.0) set the minimum RPK value to consider a marker present\n\n"
            
            "*  The list '-t clade_profiles' analysis type reports the same information of '-t marker_ab_table'\n"
            "   but the markers are reported on a clade-by-clade basis.\n"
            "$ metaphlan2.py -t clade_profiles metagenome_outfmt.bz2 --input_type bowtie2out > marker_abundance_table.txt\n\n"
            
            "*  Finally, to obtain all markers present for a specific clade and all its subclades, the \n"
            "   '-t clade_specific_strain_tracker' should be used. For example, the following command\n"
            "   is reporting the presence/absence of the markers for the B. fragulis species and its strains\n"
            "   the optional argument --min_ab specifies the minimum clade abundance for reporting the markers\n\n"
            "$ metaphlan2.py -t clade_specific_strain_tracker --clade s__Bacteroides_fragilis metagenome_outfmt.bz2 --input_type bowtie2out > marker_abundance_table.txt\n"
            
            "\n------------------------------------------------------------------- \n\n"
            "",
            formatter_class=ap.RawTextHelpFormatter,
            add_help=False )
    arg = p.add_argument

    arg( 'inp', metavar='INPUT_FILE', type=str, nargs='?', default=None, help= 
         "the input file can be:\n"
         "* a fastq file containing metagenomic reads\n"
         "OR\n"
         "* a BowTie2 produced SAM file. \n"
         "OR\n"
         "* an intermediary mapping file of the metagenome generated by a previous MetaPhlAn run \n"
         "If the input file is missing, the script assumes that the input is provided using the standard \n"
         "input, or named pipes.\n"
         "IMPORTANT: the type of input needs to be specified with --input_type" )   
    
    arg( 'output', metavar='OUTPUT_FILE', type=str, nargs='?', default=None,
         help= "the tab-separated output file of the predicted taxon relative abundances \n"
               "[stdout if not present]")


    g = p.add_argument_group('Required arguments')
    arg = g.add_argument
    input_type_choices = ['fastq','fasta','multifasta','multifastq','bowtie2out','sam'] # !!!!
    arg( '--input_type', choices=input_type_choices, required = 'True', help =  
         "set whether the input is the multifasta file of metagenomic reads or \n"
         "the SAM file of the mapping of the reads against the MetaPhlAn db.\n"
         "[default 'automatic', i.e. the script will try to guess the input format]\n" )
   
    g = p.add_argument_group('Mapping arguments')
    arg = g.add_argument
    arg( '--mpa_pkl', type=str,
         default=os.path.join(metaphlan2_script_install_folder,"db_v20","mpa_v20_m200.pkl"), 
         help = "the metadata pickled MetaPhlAn file")
    arg( '--bowtie2db', metavar="METAPHLAN_BOWTIE2_DB", type=str,
         default = os.path.join(metaphlan2_script_install_folder,"db_v20","mpa_v20_m200"),
         help = "The BowTie2 database file of the MetaPhlAn database. \n"
                "Used if --input_type is fastq, fasta, multifasta, or multifastq")
    bt2ps = ['sensitive','very-sensitive','sensitive-local','very-sensitive-local']
    arg( '--bt2_ps', metavar="BowTie2 presets", default='very-sensitive', choices=bt2ps,
         help = "presets options for BowTie2 (applied only when a multifasta file is provided)\n"
                "The choices enabled in MetaPhlAn are:\n"
                " * sensitive\n"
                " * very-sensitive\n"
                " * sensitive-local\n"
                " * very-sensitive-local\n"
                "[default very-sensitive]\n"   )
    arg( '--bowtie2_exe', type=str, default = None, help =
         'Full path and name of the BowTie2 executable. This option allows \n'
         'MetaPhlAn to reach the executable even when it is not in the system \n'
         'PATH or the system PATH is unreachable\n' )
    arg( '--bowtie2out', metavar="FILE_NAME", type=str, default = None, help = 
         "The file for saving the output of BowTie2\n" )
    arg( '--no_map', action='store_true', help=
         "Avoid storing the --bowtie2out map file\n" )
    arg( '--tmp_dir', metavar="", default=None, type=str, help = 
         "the folder used to store temporary files \n"
         "[default is the OS dependent tmp dir]\n"   )
    
    
    g = p.add_argument_group('Post-mapping arguments')
    arg = g.add_argument
    stat_choices = ['avg_g','avg_l','tavg_g','tavg_l','wavg_g','wavg_l','med']
    arg( '--tax_lev', metavar='TAXONOMIC_LEVEL', type=str, 
         choices='a'+tax_units, default='a', help = 
         "The taxonomic level for the relative abundance output:\n"
         "'a' : all taxonomic levels\n"
         "'k' : kingdoms\n"
         "'p' : phyla only\n"
         "'c' : classes only\n"
         "'o' : orders only\n"
         "'f' : families only\n"
         "'g' : genera only\n"
         "'s' : species only\n"
         "[default 'a']" )
    arg( '--min_cu_len', metavar="", default="2000", type=int, help =
         "minimum total nucleotide length for the markers in a clade for\n"
         "estimating the abundance without considering sub-clade abundances\n"
         "[default 2000]\n"   )
    arg( '--min_alignment_len', metavar="", default=None, type=int, help =
         "The sam records for aligned reads with the longest subalignment\n"
         "length smaller than this threshold will be discarded.\n"
         "[default None]\n"   )
    arg( '--ignore_viruses', action='store_true', help=
         "Do not profile viral organisms" )
    arg( '--ignore_eukaryotes', action='store_true', help=
         "Do not profile eukaryotic organisms" )
    arg( '--ignore_bacteria', action='store_true', help=
         "Do not profile bacterial organisms" )
    arg( '--ignore_archaea', action='store_true', help=
         "Do not profile archeal organisms" )
    arg( '--stat_q', metavar="", type = float, default=0.1, help = 
         "Quantile value for the robust average\n"
         "[default 0.1]"   )
    arg( '--ignore_markers', type=str, default = None, help = 
         "File containing a list of markers to ignore. \n")
    arg( '--avoid_disqm', action="store_true", help = 
         "Deactivate the procedure of disambiguating the quasi-markers based on the \n"
         "marker abundance pattern found in the sample. It is generally recommended \n"
         "too keep the disambiguation procedure in order to minimize false positives\n")
    arg( '--stat', metavar="", choices=stat_choices, default="tavg_g", type=str, help = 
         "EXPERIMENTAL! Statistical approach for converting marker abundances into clade abundances\n"
         "'avg_g'  : clade global (i.e. normalizing all markers together) average\n"
         "'avg_l'  : average of length-normalized marker counts\n"
         "'tavg_g' : truncated clade global average at --stat_q quantile\n"
         "'tavg_l' : trunated average of length-normalized marker counts (at --stat_q)\n"
         "'wavg_g' : winsorized clade global average (at --stat_q)\n"
         "'wavg_l' : winsorized average of length-normalized marker counts (at --stat_q)\n"
         "'med'    : median of length-normalized marker counts\n"
         "[default tavg_g]"   ) 
    
    arg = p.add_argument


    
    g = p.add_argument_group('Additional analysis types and arguments')
    arg = g.add_argument
    analysis_types = ['rel_ab', 'rel_ab_w_read_stats', 'reads_map', 'clade_profiles', 'marker_ab_table', 'marker_counts', 'marker_pres_table', 'clade_specific_strain_tracker']
    arg( '-t', metavar='ANALYSIS TYPE', type=str, choices = analysis_types, 
         default='rel_ab', help = 
         "Type of analysis to perform: \n"
         " * rel_ab: profiling a metagenomes in terms of relative abundances\n"
         " * rel_ab_w_read_stats: profiling a metagenomes in terms of relative abundances and estimate the number of reads comming from each clade.\n"
         " * reads_map: mapping from reads to clades (only reads hitting a marker)\n"
         " * clade_profiles: normalized marker counts for clades with at least a non-null marker\n"
         " * marker_ab_table: normalized marker counts (only when > 0.0 and normalized by metagenome size if --nreads is specified)\n"
         " * marker_counts: non-normalized marker counts [use with extreme caution]\n"
         " * marker_pres_table: list of markers present in the sample (threshold at 1.0 if not differently specified with --pres_th\n"
         "[default 'rel_ab']" )
    arg( '--nreads', metavar="NUMBER_OF_READS", type=int, default = None, help =
         "The total number of reads in the original metagenome. It is used only when \n"
         "-t marker_table is specified for normalizing the length-normalized counts \n"
         "with the metagenome size as well. No normalization applied if --nreads is not \n"
         "specified" )
    arg( '--pres_th', metavar="PRESENCE_THRESHOLD", type=int, default = 1.0, help =
         'Threshold for calling a marker present by the -t marker_pres_table option' )
    arg( '--clade', metavar="", default=None, type=str, help = 
         "The clade for clade_specific_strain_tracker analysis\n"  )
    arg( '--min_ab', metavar="", default=0.1, type=float, help = 
         "The minimum percentage abundace for the clade in the clade_specific_strain_tracker analysis\n"  )
    arg( "-h", "--help", action="help", help="show this help message and exit")

    g = p.add_argument_group('Output arguments')
    arg = g.add_argument
    arg( '-o', '--output_file',  metavar="output file", type=str, default=None, help = 
         "The output file (if not specified as positional argument)\n")
    arg('--sample_id_key',  metavar="name", type=str, default="#SampleID", 
        help =("Specify the sample ID key for this analysis."
               " Defaults to '#SampleID'."))
    arg('--sample_id',  metavar="value", type=str, 
        default="Metaphlan2_Analysis",
        help =("Specify the sample ID for this analysis."
               " Defaults to 'Metaphlan2_Analysis'."))
    arg( '-s', '--samout', metavar="sam_output_file",
        type=str, default=None, help="The sam output file\n")
    #*************************************************************
    #* Parameters related to biom file generation                *
    #*************************************************************         
    arg( '--biom', '--biom_output_file',  metavar="biom_output", type=str, default=None, help = 
         "If requesting biom file output: The name of the output file in biom format \n")

    arg( '--mdelim', '--metadata_delimiter_char',  metavar="mdelim", type=str, default="|", help = 
         "Delimiter for bug metadata: - defaults to pipe. e.g. the pipe in k__Bacteria|p__Proteobacteria \n")
    #*************************************************************
    #* End parameters related to biom file generation            *
    #*************************************************************    
    
    g = p.add_argument_group('Other arguments')
    arg = g.add_argument
    arg( '--nproc', metavar="N", type=int, default=1, help = 
         "The number of CPUs to use for parallelizing the mapping\n"
         "[default 1, i.e. no parallelism]\n" ) 
    arg( '-v','--version', action='version', version="MetaPhlAn version "+__version__+"\t("+__date__+")",
         help="Prints the current MetaPhlAn version and exit\n" )
    

    return vars(p.parse_args()) 

def run_bowtie2(  fna_in, outfmt6_out, bowtie2_db, preset, nproc, 
                  file_format = "multifasta", exe = None, 
                  samout = None,
                  min_alignment_len = None,
                  ):
    try:
        if not fna_in: # or stat.S_ISFIFO(os.stat(fna_in).st_mode):
            fna_in = "-"
        bowtie2_cmd = [ exe if exe else 'bowtie2', 
                        "--quiet", "--no-unal", 
                        "--"+preset,
                        "-S","-",
                        "-x", bowtie2_db,
                         ] + ([] if int(nproc) < 2 else ["-p",str(nproc)])
        bowtie2_cmd += ["-U", fna_in] # if not stat.S_ISFIFO(os.stat(fna_in).st_mode) else []
        bowtie2_cmd += (["-f"] if file_format == "multifasta" else []) 
        p = subp.Popen( bowtie2_cmd, stdout=subp.PIPE ) 
        lmybytes, outf = (mybytes,bz2.BZ2File(outfmt6_out, "w")) if outfmt6_out.endswith(".bz2") else (str,open( outfmt6_out, "w" ))
        
        try:
            if samout:
                if samout[-4:] == '.bz2':
                    sam_file = bz2.BZ2File(samout, 'w')
                else:
                    sam_file = open(samout, 'w')
        except IOError:
            sys.stderr.write( "IOError: Unable to open sam output file.\n" )
            sys.exit(1)

        for line in p.stdout:
            if samout:
                sam_file.write(line)
            if line[0] != '@':
                o = read_and_split_line(line)
                if o[2][-1] != '*':
                    if min_alignment_len == None\
                        or max([int(x.strip('M')) for x in\
                                re.findall(r'(\d*M)', o[5])]) >= min_alignment_len:
                        outf.write( lmybytes("\t".join([o[0],o[2]]) +"\n") )
        #if  float(sys.version_info[0]) >= 3: 
        #    for o in read_and_split(p.stdout):
        #        if o[2][-1] != '*':
        #            outf.write( bytes("\t".join([o[0],o[2]]) +"\n",encoding='utf-8') )
        #else:
        #    for o in read_and_split(p.stdout):
        #        if o[2][-1] != '*':
        #            outf.write( "\t".join([o[0],o[2]]) +"\n" )
        outf.close()
        if samout:
            sam_file.close()
        p.wait()


    except OSError:
        sys.stderr.write( "OSError: fatal error running BowTie2. Is BowTie2 in the system path?\n" )
        sys.exit(1)
    except ValueError:
        sys.stderr.write( "ValueError: fatal error running BowTie2.\n" )
        sys.exit(1)
    except IOError:
        sys.stderr.write( "IOError: fatal error running BowTie2.\n" )
        sys.exit(1)
    if p.returncode == 13:
        sys.stderr.write( "Permission Denied Error: fatal error running BowTie2." 
          "Is the BowTie2 file in the path with execution and read permissions?\n" )
        sys.exit(1)
    elif p.returncode != 0:
        sys.stderr.write("Error while running bowtie2.\n")
        sys.exit(1)

#def guess_input_format( inp_file ):
#    if "," in inp_file:
#        sys.stderr.write( "Sorry, I cannot guess the format of the input, when "
#                          "more than one file is specified. Please set the --input_type parameter \n" )
#        sys.exit(1) 
#
#    with open( inp_file ) as inpf:
#        for i,l in enumerate(inpf):
#            line = l.strip()
#            if line[0] == '#': continue
#            if line[0] == '>': return 'multifasta'
#            if line[0] == '@': return 'multifastq'
#            if len(l.split('\t')) == 2: return 'bowtie2out'
#            if i > 20: break
#    return None

class TaxClade:
    min_cu_len = -1
    markers2lens = None
    stat = None
    quantile = None
    avoid_disqm = False

    def __init__( self, name, uncl = False, id_int = 0 ):
        self.children, self.markers2nreads = {}, {}
        self.name, self.father = name, None
        self.uncl, self.subcl_uncl = uncl, False
        self.abundance, self.uncl_abundance = None, 0 
        self.id = id_int

    def add_child( self, name, id_int ):
        new_clade = TaxClade( name, id_int=id_int )
        self.children[name] = new_clade
        new_clade.father = self
        return new_clade

    
    def get_terminals( self ):
        terms = []
        if not self.children:
            return [self]
        for c in self.children.values():
            terms += c.get_terminals()
        return terms


    def get_full_name( self ):
        fullname = [self.name]
        cl = self.father
        while cl:
            fullname = [cl.name] + fullname
            cl = cl.father
        return "|".join(fullname[1:])

    def get_normalized_counts( self ):
        return [(m,float(n)*1000.0/self.markers2lens[m]) 
                    for m,n in self.markers2nreads.items()]

    def compute_abundance( self ):
        if self.abundance is not None: return self.abundance
        sum_ab = sum([c.compute_abundance() for c in self.children.values()]) 
        rat_nreads = sorted([(self.markers2lens[m],n) 
                                    for m,n in self.markers2nreads.items()],
                                            key = lambda x: x[1])

        rat_nreads, removed = [], []
        for m,n in self.markers2nreads.items():
            misidentified = False

            if not self.avoid_disqm:
                for e in self.markers2exts[m]:
                    toclade = self.taxa2clades[e]
                    m2nr = toclade.markers2nreads
                    tocladetmp = toclade
                    while len(tocladetmp.children) == 1:
                        tocladetmp = list(tocladetmp.children.values())[0]
                        m2nr = tocladetmp.markers2nreads
    
                    nonzeros = sum([v>0 for v in m2nr.values()])
                    if len(m2nr):
                        if float(nonzeros) / len(m2nr) > 0.33:
                            misidentified = True
                            removed.append( (self.markers2lens[m],n) )
                            break
            if not misidentified:
                rat_nreads.append( (self.markers2lens[m],n) ) 
       
        if not self.avoid_disqm and len(removed):
            n_rat_nreads = float(len(rat_nreads))
            n_removed = float(len(removed))
            n_tot = n_rat_nreads + n_removed
            n_ripr = 10
            
            if len(self.get_terminals()) < 2:
                n_ripr = 0

            if "k__Viruses" in self.get_full_name():
                n_ripr = 0

            if n_rat_nreads < n_ripr and n_tot > n_rat_nreads:
                rat_nreads += removed[:n_ripr-int(n_rat_nreads)]

        
        rat_nreads = sorted(rat_nreads, key = lambda x: x[1])

        rat_v,nreads_v = zip(*rat_nreads) if rat_nreads else ([],[])
        rat, nrawreads, loc_ab = float(sum(rat_v)) or -1.0, sum(nreads_v), 0.0
        quant = int(self.quantile*len(rat_nreads))
        ql,qr,qn = (quant,-quant,quant) if quant else (None,None,0)
     
        if self.name[0] == 't' and (len(self.father.children) > 1 or "_sp" in self.father.name or "k__Viruses" in self.get_full_name()):
            non_zeros = float(len([n for r,n in rat_nreads if n > 0])) 
            nreads = float(len(rat_nreads))
            if nreads == 0.0 or non_zeros / nreads < 0.7:
                self.abundance = 0.0
                return 0.0

        if rat < 0.0:
            pass
        elif self.stat == 'avg_g' or (not qn and self.stat in ['wavg_g','tavg_g']):
            loc_ab = nrawreads / rat if rat >= 0 else 0.0
        elif self.stat == 'avg_l' or (not qn and self.stat in ['wavg_l','tavg_l']):
            loc_ab = np.mean([float(n)/r for r,n in rat_nreads]) 
        elif self.stat == 'tavg_g':
            wnreads = sorted([(float(n)/r,r,n) for r,n in rat_nreads], key=lambda x:x[0])
            den,num = zip(*[v[1:] for v in wnreads[ql:qr]])
            loc_ab = float(sum(num))/float(sum(den)) if any(den) else 0.0
        elif self.stat == 'tavg_l':
            loc_ab = np.mean(sorted([float(n)/r for r,n in rat_nreads])[ql:qr])
        elif self.stat == 'wavg_g':
            vmin, vmax = nreads_v[ql], nreads_v[qr]
            wnreads = [vmin]*qn+list(nreads_v[ql:qr])+[vmax]*qn
            loc_ab = float(sum(wnreads)) / rat  
        elif self.stat == 'wavg_l':
            wnreads = sorted([float(n)/r for r,n in rat_nreads])
            vmin, vmax = wnreads[ql], wnreads[qr]
            wnreads = [vmin]*qn+list(wnreads[ql:qr])+[vmax]*qn
            loc_ab = np.mean(wnreads) 
        elif self.stat == 'med':
            loc_ab = np.median(sorted([float(n)/r for r,n in rat_nreads])[ql:qr]) 
        
        self.abundance = loc_ab
        if rat < self.min_cu_len and self.children:
            self.abundance = sum_ab
        elif loc_ab < sum_ab:
            self.abundance = sum_ab

        if self.abundance > sum_ab and self.children: # *1.1??
            self.uncl_abundance = self.abundance - sum_ab
        self.subcl_uncl = not self.children and self.name[0] not in tax_units[-2:] 

        return self.abundance

    def get_all_abundances( self ):
        ret = [(self.name,self.abundance)]
        if self.uncl_abundance > 0.0:
            lchild = list(self.children.values())[0].name[:3]
            ret += [(lchild+self.name[3:]+"_unclassified",self.uncl_abundance)]
        if self.subcl_uncl and self.name[0] != tax_units[-2]:
            cind = tax_units.index( self.name[0] )
            ret += [(   tax_units[cind+1]+self.name[1:]+"_unclassified",
                        self.abundance)]
        for c in self.children.values():
            ret += c.get_all_abundances()
        return ret


class TaxTree:
    def __init__( self, mpa, markers_to_ignore = None ): #, min_cu_len ):
        self.root = TaxClade( "root" )
        self.all_clades, self.markers2lens, self.markers2clades, self.taxa2clades, self.markers2exts = {}, {}, {}, {}, {}
        TaxClade.markers2lens = self.markers2lens
        TaxClade.markers2exts = self.markers2exts
        TaxClade.taxa2clades = self.taxa2clades
        self.id_gen = itertools.count(1)

        clades_txt = ((l.strip().split("|"),n) for l,n in mpa_pkl['taxonomy'].items())        
        for clade,lenc in clades_txt:
            father = self.root
            for clade_lev in clade: # !!!!! [:-1]:
                if not clade_lev in father.children:
                    father.add_child( clade_lev, id_int=next(self.id_gen) )
                    self.all_clades[clade_lev] = father.children[clade_lev]
                if clade_lev[0] == "t":
                    self.taxa2clades[clade_lev[3:]] = father

                father = father.children[clade_lev]
                if clade_lev[0] == "t":
                    father.glen = lenc

        def add_lens( node ):
            if not node.children:
                return node.glen
            lens = [] 
            for c in node.children.values():
                lens.append( add_lens( c ) )
            node.glen = sum(lens) / len(lens)
            return node.glen
        add_lens( self.root )
        
        for k,p in mpa_pkl['markers'].items():
            if k in markers_to_exclude:
                continue
            if k in markers_to_ignore:
                continue
            self.markers2lens[k] = p['len']
            self.markers2clades[k] = p['clade']
            self.add_reads( k, 0  )
            self.markers2exts[k] = p['ext']

    def set_min_cu_len( self, min_cu_len ):
        TaxClade.min_cu_len = min_cu_len

    def set_stat( self, stat, quantile, avoid_disqm = False ):
        TaxClade.stat = stat
        TaxClade.quantile = quantile
        TaxClade.avoid_disqm = avoid_disqm

    def add_reads(  self, marker, n, 
                    ignore_viruses = False, ignore_eukaryotes = False, 
                    ignore_bacteria = False, ignore_archaea = False  ):
        clade = self.markers2clades[marker]
        cl = self.all_clades[clade]
        if ignore_viruses or ignore_eukaryotes or ignore_bacteria or ignore_archaea:
            cn = cl.get_full_name()
            if ignore_viruses and cn.startswith("k__Viruses"):
                return ""
            if ignore_eukaryotes and cn.startswith("k__Eukaryota"):
                return ""
            if ignore_archaea and cn.startswith("k__Archaea"):
                return ""
            if ignore_bacteria and cn.startswith("k__Bacteria"):
                return ""
        while len(cl.children) == 1:
            cl = list(cl.children.values())[0]
        cl.markers2nreads[marker] = n
        return cl.get_full_name()
   

    def markers2counts( self ):
        m2c = {}
        for k,v in self.all_clades.items():
            for m,c in v.markers2nreads.items():
                m2c[m] = c
        return m2c

    def clade_profiles( self, tax_lev, get_all = False  ):
        cl2pr = {}
        for k,v in self.all_clades.items():
            if tax_lev and not k.startswith(tax_lev): 
                continue
            prof = v.get_normalized_counts()
            if not get_all and ( len(prof) < 1 or not sum([p[1] for p in prof]) > 0.0 ):
                continue
            cl2pr[v.get_full_name()] = prof
        return cl2pr
            
    def relative_abundances( self, tax_lev  ):
        cl2ab_n = dict([(k,v) for k,v in self.all_clades.items() 
                    if k.startswith("k__") and not v.uncl])
     
        cl2ab, cl2glen, tot_ab = {}, {}, 0.0 
        for k,v in cl2ab_n.items():
            tot_ab += v.compute_abundance()

        for k,v in cl2ab_n.items():
            for cl,ab in v.get_all_abundances():
                if not tax_lev:
                    if cl not in self.all_clades:
                        to = tax_units.index(cl[0])
                        t = tax_units[to-1]
                        cl = t + cl.split("_unclassified")[0][1:]
                        cl = self.all_clades[cl].get_full_name()
                        spl = cl.split("|")
                        cl = "|".join(spl+[tax_units[to]+spl[-1][1:]+"_unclassified"])
                        glen = self.all_clades[spl[-1]].glen
                    else:
                        glen = self.all_clades[cl].glen
                        cl = self.all_clades[cl].get_full_name() 
                elif not cl.startswith(tax_lev):
                    if cl in self.all_clades:
                        glen = self.all_clades[cl].glen
                    else:
                        glen = 1.0
                    continue
                cl2ab[cl] = ab
                cl2glen[cl] = glen 

        ret_d = dict([( k, float(v) / tot_ab if tot_ab else 0.0) for k,v in cl2ab.items()])
        ret_r = dict([( k, (v,cl2glen[k],float(v)*cl2glen[k])) for k,v in cl2ab.items()])
        #ret_r = dict([( k, float(v) / tot_ab if tot_ab else 0.0) for k,v in cl2ab.items()])
        if tax_lev:
            ret_d[tax_lev+"unclassified"] = 1.0 - sum(ret_d.values())
        return ret_d, ret_r

def map2bbh( mapping_f, input_type = 'bowtie2out', min_alignment_len = None):
    if not mapping_f:
        ras, ras_line, inpf = plain_read_and_split, plain_read_and_split_line, sys.stdin
    else:
        if mapping_f.endswith(".bz2"):
            ras, ras_line, inpf = read_and_split, read_and_split_line, bz2.BZ2File( mapping_f, "r" )
        else:
            ras, ras_line, inpf = plain_read_and_split,\
                                  plain_read_and_split_line,\
                                  open( mapping_f )

    reads2markers, reads2maxb = {}, {}
    if input_type == 'bowtie2out':
        for r,c in ras(inpf):
            reads2markers[r] = c
    elif input_type == 'sam':
        for line in inpf:
            o = ras_line(line)
            if o[0][0] != '@' and o[2][-1] != '*':
                if min_alignment_len == None\
                    or max([int(x.strip('M')) for x in\
                            re.findall(r'(\d*M)', o[5])]) >= min_alignment_len:
                    reads2markers[o[0]] = o[2]
    inpf.close()

    markers2reads = defdict( set )
    for r,m in reads2markers.items():
        markers2reads[m].add( r )

    return markers2reads
    
    
def maybe_generate_biom_file(pars, abundance_predictions):
    if not pars['biom']:
        return None
    if not abundance_predictions:
        return open(pars['biom'], 'w').close()

    delimiter = "|" if len(pars['mdelim']) > 1 else pars['mdelim']
    def istip(clade_name):
        end_name = clade_name.split(delimiter)[-1]
        return end_name.startswith("t__") or end_name.endswith("_unclassified")

    def findclade(clade_name):
        if clade_name.endswith('_unclassified'):
            name = clade_name.split(delimiter)[-2]
        else:
            name = clade_name.split(delimiter)[-1]
        return tree.all_clades[name]

    def to_biomformat(clade_name):
        return { 'taxonomy': clade_name.split(delimiter) }

    clades = iter( (abundance, findclade(name)) 
                   for (name, abundance) in abundance_predictions
                   if istip(name) )
    packed = iter( ([abundance], clade.get_full_name(), clade.id)
                   for (abundance, clade) in clades )

    #unpack that tuple here to stay under 80 chars on a line
    data, clade_names, clade_ids = zip(*packed)
    # biom likes column vectors, so we give it an array like this:
    # np.array([a],[b],[c])
    data = np.array(data)
    sample_ids = [pars['sample_id']]
    table_id='MetaPhlAn2_Analysis'
    json_key = "MetaPhlAn2"
  


    #**********************************************************************************************
    #  Modification of Code :                                                                     *
    #  Modified the code so instead of using the current clade IDs, which are numbers, we will    *
    #      use the clade_names                                                                    *
    #      Users reported the biom output is invalid and also the IDs were changing from run to   *
    #      run.                                                                                   *
    #  George Weingart    05/22/2017   george.weingart@mail.com                                   *
    #**********************************************************************************************
    if LooseVersion(biom.__version__) < LooseVersion("2.0.0"):
        biom_table = biom.table.table_factory(
            data, 
	    sample_ids,
            ######## clade_ids,     #Modified by George Weingart 5/22/2017 - We will use instead the clade_names
            clade_names,            #Modified by George Weingart 5/22/2017 - We will use instead the clade_names
            sample_metadata      = None,
            observation_metadata = map(to_biomformat, clade_names),
            table_id             = table_id,
            constructor          = biom.table.DenseOTUTable
        )
        with open(pars['biom'], 'w') as outfile:
            json.dump( biom_table.getBiomFormatObject(json_key),
                           outfile )
    else:  # Below is the biom2 compatible code
        biom_table = biom.table.Table(
            data, 
            #clade_ids,           #Modified by George Weingart 5/22/2017 - We will use instead the clade_names
            clade_names,          #Modified by George Weingart 5/22/2017 - We will use instead the clade_names 
            sample_ids,
            sample_metadata      = None,
            observation_metadata = map(to_biomformat, clade_names),
            table_id             = table_id,
            input_is_dense       = True
        )
        
        with open(pars['biom'], 'w') as outfile:  
            biom_table.to_json( json_key,
                                direct_io = outfile )

    return True


if __name__ == '__main__':
    pars = read_params( sys.argv )    
    #if pars['inp'] is None and ( pars['input_type'] is None or  pars['input_type'] == 'automatic'): 
    #    sys.stderr.write( "The --input_type parameter need top be specified when the "
    #                      "input is provided from the standard input.\n"
    #                      "Type metaphlan.py -h for more info\n")
    #    sys.exit(0)

    if pars['bt2_ps'] in [
                          "sensitive-local",
                          "very-sensitive-local"
                          ]\
        and pars['min_alignment_len'] == None:
            pars['min_alignment_len'] = 100
            sys.stderr.write('Warning! bt2_ps is set to local mode, '\
                             'and min_alignment_len is None, '
                             'I automatically set min_alignment_len to 100! '\
                             'If you do not like, rerun the command and set '\
                             'min_alignment_len to a specific value.\n'
                            )

    if pars['input_type'] == 'fastq':
        pars['input_type'] = 'multifastq'
    if pars['input_type'] == 'fasta':
        pars['input_type'] = 'multifasta'

    #if pars['input_type'] == 'automatic':
    #    pars['input_type'] = guess_input_format( pars['inp'] )
    #    if not pars['input_type']:
    #        sys.stderr.write( "Sorry, I cannot guess the format of the input file, please "
    #                          "specify the --input_type parameter \n" )
    #        sys.exit(1) 

    # check for the mpa_pkl file
    if not os.path.isfile(pars['mpa_pkl']):
        sys.stderr.write("Error: Unable to find the mpa_pkl file at: " + pars['mpa_pkl'] +
                         "\nExpecting location ${mpa_dir}/db_v20/map_v20_m200.pkl "
                         "\nSelect the file location with the option --mpa_pkl.\n"
                         "Exiting...\n\n")
        sys.exit(1)           

    if pars['ignore_markers']:
        with open(pars['ignore_markers']) as ignv:
            ignore_markers = set([l.strip() for l in ignv])
    else:
        ignore_markers = set()

    no_map = False
    if pars['input_type'] == 'multifasta' or pars['input_type'] == 'multifastq':
        bow = pars['bowtie2db'] is not None
        if not bow:
            sys.stderr.write( "No MetaPhlAn BowTie2 database provided\n "
                              "[--bowtie2db options]!\n"
                              "Exiting...\n\n" )
            sys.exit(1)
        if pars['no_map']:
            pars['bowtie2out'] = tf.NamedTemporaryFile(dir=pars['tmp_dir']).name
            no_map = True
        else:
            if bow and not pars['bowtie2out']:
                if pars['inp'] and "," in  pars['inp']:
                    sys.stderr.write( "Error! --bowtie2out needs to be specified when multiple "
                                      "fastq or fasta files (comma separated) are provided"  )
                    sys.exit(1)
                fname = pars['inp']
                if fname is None:
                    fname = "stdin_map"
                elif stat.S_ISFIFO(os.stat(fname).st_mode):
                    fname = "fifo_map"
                pars['bowtie2out'] = fname + ".bowtie2out.txt"

            if os.path.exists( pars['bowtie2out'] ):
                sys.stderr.write(   
                    "BowTie2 output file detected: " + pars['bowtie2out'] + "\n"
                    "Please use it as input or remove it if you want to "
                    "re-perform the BowTie2 run.\n"
                    "Exiting...\n\n" )
                sys.exit(1)

        if bow and not all([os.path.exists(".".join([str(pars['bowtie2db']),p]))
                        for p in ["1.bt2", "2.bt2", "3.bt2","4.bt2","1.bt2","2.bt2"]]):
            sys.stderr.write( "No MetaPhlAn BowTie2 database found "
                              "[--bowtie2db option]! "
                              "(or wrong path provided)."
                              "\nExpecting location ${mpa_dir}/db_v20/map_v20_m200 "
                              "\nExiting... " )
            sys.exit(1)

        if bow:
            run_bowtie2( pars['inp'], pars['bowtie2out'], pars['bowtie2db'], 
                         pars['bt2_ps'], pars['nproc'], file_format = pars['input_type'],
                         exe = pars['bowtie2_exe'],
                         samout = pars['samout'],
                         min_alignment_len = pars['min_alignment_len'])
            pars['input_type'] = 'bowtie2out'
        
        pars['inp'] = pars['bowtie2out'] # !!!

    with open( pars['mpa_pkl'], 'rb' ) as a:
        mpa_pkl = pickle.loads( bz2.decompress( a.read() ) )

    tree = TaxTree( mpa_pkl, ignore_markers )
    tree.set_min_cu_len( pars['min_cu_len'] )
    tree.set_stat( pars['stat'], pars['stat_q'], pars['avoid_disqm']  )

    markers2reads = map2bbh( 
                            pars['inp'], 
                            pars['input_type'],
                            pars['min_alignment_len']
                            )
    if no_map:
        os.remove( pars['inp'] )         

    map_out = []
    for marker,reads in markers2reads.items():
        if marker not in tree.markers2lens:
            continue
        tax_seq = tree.add_reads( marker, len(reads), 
                                  ignore_viruses = pars['ignore_viruses'],
                                  ignore_eukaryotes = pars['ignore_eukaryotes'],
                                  ignore_bacteria = pars['ignore_bacteria'],
                                  ignore_archaea = pars['ignore_archaea'],
                                  )
        if tax_seq:
            map_out +=["\t".join([r,tax_seq]) for r in reads]
    
    if pars['output'] is None and pars['output_file'] is not None:
        pars['output'] = pars['output_file']

    with (open(pars['output'],"w") if pars['output'] else sys.stdout) as outf:
        outf.write('\t'.join((pars["sample_id_key"], pars["sample_id"])) + '\n')
        if pars['t'] == 'reads_map':
            outf.write( "\n".join( map_out ) + "\n" )
        elif pars['t'] == 'rel_ab':
            cl2ab, _ = tree.relative_abundances( 
                        pars['tax_lev']+"__" if pars['tax_lev'] != 'a' else None )
            outpred = [(k,round(v*100.0,5)) for k,v in cl2ab.items() if v > 0.0]
            if outpred:
                for k,v in sorted(  outpred, reverse=True,
                                    key=lambda x:x[1]+(100.0*(8-x[0].count("|")))  ): 
                    outf.write( "\t".join( [k,str(v)] ) + "\n" )   
            else:
                outf.write( "unclassified\t100.0\n" )
            maybe_generate_biom_file(pars, outpred)
        elif pars['t'] == 'rel_ab_w_read_stats':
            cl2ab, rr = tree.relative_abundances( 
                        pars['tax_lev']+"__" if pars['tax_lev'] != 'a' else None )
            outpred = [(k,round(v*100.0,5)) for k,v in cl2ab.items() if v > 0.0]
            totl = 0
            if outpred:
                outf.write( "\t".join( [    "#clade_name",
                                            "relative_abundance",
                                            "coverage",
                                            "average_genome_length_in_the_clade",
                                            "estimated_number_of_reads_from_the_clade" ]) +"\n" )

                for k,v in sorted(  outpred, reverse=True,
                                    key=lambda x:x[1]+(100.0*(8-x[0].count("|")))  ): 
                    outf.write( "\t".join( [    k,
                                                str(v),
                                                str(rr[k][0]) if k in rr else "-",
                                                str(rr[k][1]) if k in rr else "-",
                                                str(int(round(rr[k][2],0)) if k in rr else "-")   
                                                ] ) + "\n" )   
                    if "|" not in k:
                        totl += (int(round(rr[k][2],0)) if k in rr else 0)

                outf.write( "#estimated total number of reads from known clades: " + str(totl)+"\n")
            else:
                outf.write( "unclassified\t100.0\n" )
            maybe_generate_biom_file(pars, outpred)

        elif pars['t'] == 'clade_profiles':
            cl2pr = tree.clade_profiles( pars['tax_lev']+"__" if pars['tax_lev'] != 'a' else None  )
            for c,p in cl2pr.items():
                mn,n = zip(*p)
                outf.write( "\t".join( [""]+[str(s) for s in mn] ) + "\n" )
                outf.write( "\t".join( [c]+[str(s) for s in n] ) + "\n" )
        elif pars['t'] == 'marker_ab_table':
            cl2pr = tree.clade_profiles( pars['tax_lev']+"__" if pars['tax_lev'] != 'a' else None  )
            for v in cl2pr.values():
                outf.write( "\n".join(["\t".join([str(a),str(b/float(pars['nreads'])) if pars['nreads'] else str(b)]) 
                                for a,b in v if b > 0.0]) + "\n" )
        elif pars['t'] == 'marker_pres_table':
            cl2pr = tree.clade_profiles( pars['tax_lev']+"__" if pars['tax_lev'] != 'a' else None  )
            for v in cl2pr.values():
                strout = ["\t".join([str(a),"1"]) for a,b in v if b > pars['pres_th']]
                if strout:
                    outf.write( "\n".join(strout) + "\n" )

        elif pars['t'] == 'marker_counts':
            outf.write( "\n".join( ["\t".join([m,str(c)]) for m,c in tree.markers2counts().items() ]) +"\n" )

        elif pars['t'] == 'clade_specific_strain_tracker':
            cl2pr = tree.clade_profiles( None, get_all = True  )
            cl2ab, _ = tree.relative_abundances( None )
            strout = []
            for cl,v in cl2pr.items():
                if cl.endswith(pars['clade']) and cl2ab[cl]*100.0 < pars['min_ab']:
                    strout = []
                    break
                if pars['clade'] in cl:
                    strout += ["\t".join([str(a),str(int(b > pars['pres_th']))]) for a,b in v]
            if strout:
                strout = sorted(strout,key=lambda x:x[0])
                outf.write( "\n".join(strout) + "\n" )
            else:
                sys.stderr.write("Clade "+pars['clade']+" not present at an abundance >"+str(round(pars['min_ab'],2))+"%, "
                                 "so no clade specific markers are reported\n")
