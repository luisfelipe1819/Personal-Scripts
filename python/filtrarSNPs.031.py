# coding: utf-8
from __future__ import absolute_import, division, print_function, unicode_literals

import time
import re
import argparse as ap
import math as m
from pandas import DataFrame,read_table
import sys
import os 


    
    
def metodo_matriz (table_SNPs,dic_chrm_sizePb,xyLen,cChom,cPosic,cLog,cPvl,list_chrm_sample,factor_divisor_matriz):
    inicio = time.time()

    if cLog:
        vMax = table_SNPs.iloc[:,cLog].max()
        vMin = table_SNPs.iloc[:,cLog].min()
    else:
        vMax = -m.log10(table_SNPs.iloc[:,cPvl].min())
        vMin = -m.log10(table_SNPs.iloc[:,cPvl].max())
    
    
    pbPerPixel_X = sum(dic_chrm_sizePb.values())/(xyLen[0]/factor_divisor_matriz)
    valoresPerPixel_y  = (vMax-vMin)/(xyLen[1]/factor_divisor_matriz)
    

    
    dic_chrm_pbAbsolutas = {}
    posicion_chrmAnterior = 0
    for chrm in list_chrm_sample:
        dic_chrm_pbAbsolutas[chrm] = posicion_chrmAnterior
        posicion_chrmAnterior += dic_chrm_sizePb[chrm]
    
    dic_xyCoordenadas_SNPsMayorLog = {}
    list_absPosic = []

    for num,SNP in enumerate(table_SNPs.itertuples()):

        if cLog:
            log = SNP[cLog+1]
        else:
            log = -m.log10(SNP[cPvl+1])
        chrm = str(SNP[cChom+1])
        poscAbs = dic_chrm_pbAbsolutas[chrm]+SNP[cPosic+1]
        list_absPosic.append(poscAbs)
        x = int(poscAbs/pbPerPixel_X)
        y = int((log- vMin)/valoresPerPixel_y)

        if dic_xyCoordenadas_SNPsMayorLog.get((x,y),False):
            if dic_xyCoordenadas_SNPsMayorLog[(x,y)][1] < log:
                dic_xyCoordenadas_SNPsMayorLog[(x,y)] = (SNP[0],log)
        else:
            dic_xyCoordenadas_SNPsMayorLog[(x,y)] = (SNP[0],log)


    table_SNPs['Absolute Position'] = list_absPosic
    del list_absPosic

    return (set(zip(*dic_xyCoordenadas_SNPsMayorLog.values())[0]),time.time()-inicio)
    

if __name__ == '__main__':

    inicio = time.time()
    parser = ap.ArgumentParser(description="Scrip para filtrar SNPs por -log(pValue)",
                               formatter_class=ap.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--SNPsfile', default=None, type=str, help="SNP file")
    parser.add_argument('-c', '--columChrm', default=None, type=int, help="What is the number of the column where the Chromosomes are?")
    parser.add_argument('-v', '--columPvalue', default=None, type=int, help="What is the number of the column where the Pvalues are?")
    parser.add_argument('-l', '--columLog', default=None, type=int, help="What is the number of the column where the -log(pvalues) are?")
    parser.add_argument('-p', '--columPosic', default=None, type=int, help="What is the number of the column where the SNPs positions are")
    parser.add_argument('-x', '--xLen', default=None, type=int, help="Width in pixels")
    parser.add_argument('-y', '--yLen', default=None, type=int, help="Height in pixels")
    parser.add_argument('-d', '--factorDivisor', default=4, type=float, help="Division factor to reduce the pixel canvas")

    args = parser.parse_args()

    #args.SNPsfile
    
    #args.columChrm
    #args.columPvalue
    #args.columLog
    #args.columPosic
    
    #args.xLen
    #args.yLen
    
    #args.factorDivisor
    
    #print('Peak memory = ',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
    
    xLen = args.xLen
    yLen = args.yLen
    if args.columChrm  is None or args.SNPsfile is None or args.columPosic is None or xLen is None or yLen is None:
        parser.error("Error. Need specify arguments")
    if args.columLog is None and args.columPvalue is not None:
        cLog = args.columLog
        cPvl = args.columPvalue - 1
    elif args.columPvalue is None and args.columLog is not None:
        cLog = args.columLog - 1
        cPvl = args.columPvalue
    elif args.columPvalue is not None and args.columLog is not None:
        cLog = args.columLog - 1
        cPvl = args.columPvalue - 1
    else:
        parser.error("Error. Need specify p Value or -Log column")
    
    cChom = args.columChrm - 1
    cPosic = args.columPosic - 1
    
    table_SNPs = read_table(args.SNPsfile)
    

    
    set_chrms = set([str(i) for i in table_SNPs.iloc[:,cChom].unique()])

    dic_chrm_sizePb = {'1':248956422,'2':242193529,'3':198295559,'4':190214555,'5':181538259,'6':170805979,'7':159345973,
                     '8':145138636,'9':138394717,'10':133797422,'11':135086622,'12':133275309,'13':114364328,'14':107043718,
                     '15':101991189,'16':90338345,'17':83257441,'18':80373285,'19':58617616,'20':64444167,'21':46709983,
                     '22':50818468,'X':156040895,'Y':57227415}
                     
    dic_chrm_sizePb = {chrm:dic_chrm_sizePb\
    [re.search('(\d+)',chrm).group(1) if re.search('(\d+)',chrm) else chrm.upper()] for chrm in set_chrms}
    list_chrm_sample = sorted(set_chrms,key=lambda x : int(re.search('(\d+)',x).group(1)) if re.search('(\d+)',x) else x )
    
    set_indexSNPs_aGraficar,tiempo_Funcion = metodo_matriz(table_SNPs,dic_chrm_sizePb,(xLen,yLen),cChom,cPosic,cLog,cPvl,list_chrm_sample,args.factorDivisor)


    table_SNPs.loc[set_indexSNPs_aGraficar].to_csv('{}.fltd'.format(args.SNPsfile),sep=str('\t'),index=False)
    

    print('Total time = ',time.time()-inicio)
    print('SNPs = ',len(set_indexSNPs_aGraficar),' Algorithm time = ',tiempo_Funcion)
        
        
