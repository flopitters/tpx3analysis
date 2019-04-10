#! /usr/bin/env python

# ============================================================================
# File: acq_conversion.py
# -----------------------
#
# Description:
# Clustering and sorting of Timepix3 data.
#
# Notes:
# - Data files are expected to be recorded with the SPIDR python scripts
# - Data files are expected to be named 'Fe_%s_%s_data.root' % (assembly, condition), e.g. 'Fe_W0005_E02_thr_1190_ik_10_data.root'
# 
#
# ============================================================================

from ROOT import *
from array import array
from optparse import OptionParser
from scipy.cluster.hierarchy import fclusterdata
import sys
import time
import os
import numpy as np



parser = OptionParser()
parser.add_option('-b', '--assembly', help='Assembly name', dest='ASSEMBLY')
parser.add_option('-s', '--source', help='Source name', dest='SOURCE')
parser.add_option('-c', '--condition', help='Condition name', dest='CONDITION')

(options, args) = parser.parse_args()

if(options.ASSEMBLY):
    assembly=options.ASSEMBLY

else:
    print "Please specify assembly"
    parser.print_help()
    exit()

if(options.SOURCE):
    source=options.SOURCE

else:
    print "Please specify source"
    parser.print_help()
    exit()

if(options.CONDITION):
    cond=options.CONDITION

else:
    cond = ""
#    print "Please specify conditions"
#    parser.print_help()
#    exit()

  
last_time = time.time()
all_clusters = []


def ReadFile(assembly, source):
    home_path = os.environ['HOME']
    data_path = "%s/Documents/Data/tpx3/%s/%s" % (home_path, assembly, source)
    
    assembly_start = assembly.split("-")[0]

    #file_list = ['data_0.dat']
    #file_list = ['data_0.dat', 'data_1.dat', 'data_2.dat']
    file_list = ['data_0.dat', 'data_1.dat', 'data_2.dat', 'data_3.dat', 'data_4.dat', 'data_5.dat', 'data_6.dat', 'data_7.dat', 'data_8.dat']
    #file_list = ['data_0.dat', 'data_1.dat', 'data_2.dat', 'data_3.dat', 'data_4.dat', 'data_5.dat', 'data_6.dat', 'data_7.dat', 'data_8.dat', 'data_9.dat']
    if source == "fe":
        outfile = TFile("%s/Fe_%s_%s_data.root" % (data_path, assembly, cond), "RECREATE")
    elif source == "mn": # fe55 source
        outfile = TFile("%s/Mn_%s_%s_data.root" % (data_path, assembly, cond), "RECREATE")
    elif source == "am241": # am241 source
        outfile = TFile("%s/Am241_%s_%s_data.root" % (data_path, assembly, cond), "RECREATE")
    elif source == "cd":
        outfile = TFile("%s/Cd_%s_%s_data.root" % (data_path, assembly, cond), "RECREATE")
    elif source == "ag": # cd109 source
        outfile = TFile("%s/Ag_%s_%s_data.root" % (data_path, assembly, cond), "RECREATE")
    elif source == "in":   
        outfile = TFile("%s/In_%s_%s_data.root" % (data_path, assembly, cond), "RECREATE") 
    elif source == "cu":
        outfile = TFile("%s/Cu_%s_%s_data.root" % (data_path, assembly, cond), "RECREATE")
    elif source == "pb":
        outfile = TFile("%s/Pb_%s_%s_data.root" % (data_path, assembly, cond), "RECREATE")
    elif source == "unknown":
        outfile = TFile("%s/Unknown_%s_%s_data.root" % (data_path, assembly, cond), "RECREATE")
       
       
    ## Initialise             
    X = []
    Y = []
    TOT = []
    TOT_HD = []
    
    nFrames = 0
    size_cluster = 0
    sizes = []

    pixel_1hit = 0
    pixel_2hit = 0
    pixel_3hit = 0
    pixel_4hit = 0

    hist0 = TH1F("clustersize", "Clustersize", 6, 0, 6) # 
    hist1 = TH2F("hitmap", "Hit Map", 256, 0, 256, 256, 0, 256)
    hist2 = TH2F("hitmap - filtered", "Hit Map Filtered", 256, 0, 256, 256, 0, 256) 
    hist3 = TH1F("hits per frame", " Hits per Frame", 50, 0, 500) # hits per frame

    size1cnt = 0
    allcnt = 0
    everycnt = 0
    headercnt = 1
    linecnt = 0
    last_time = time.time()

    singleHits = TTree("singleHits", "single hit Clusters")
    allHits = TTree("allHits", "all hit Clusters")
    everyPix = TTree("everyPix", "every single Pixel")

    xt = array('i', [ 0 ])
    yt = array('i', [ 0 ])
    tott = array('f', [ 0 ])
    tott_hd = array('f', [ 0 ])
    csize = array('i', [ 0 ])

    singleHits.Branch('col', xt, 'col/I')
    singleHits.Branch('row', yt, 'row/I')
    singleHits.Branch('tot', tott, 'tot/F')
    singleHits.Branch('tot_hd', tott_hd, 'tot_hd/F')
    singleHits.Branch('clustersize', csize, 'size/I')

    xt_all = array('i', [ 0 ])
    yt_all = array('i', [ 0 ])
    tott_all = array('f', [ 0 ])
    tott_hd_all = array('f', [ 0 ])
    csize_all = array('i', [ 0 ])

    allHits.Branch('col', xt_all, 'col/I')
    allHits.Branch('row', yt_all, 'row/I')
    allHits.Branch('tot', tott_all, 'tot/F')
    allHits.Branch('tot_hd', tott_hd_all, 'tot_hd/F')
    allHits.Branch('clustersize', csize_all, 'size/I')

    xt_every = array('i', [ 0 ])
    yt_every = array('i', [ 0 ])
    tott_every = array('f', [ 0 ])
    tott_hd_every = array('f', [ 0 ])

    everyPix.Branch('col', xt_every, 'col/I')
    everyPix.Branch('row', yt_every, 'row/I')
    everyPix.Branch('tot', tott_every, 'tot/F')
    everyPix.Branch('tot_hd', tott_hd_every, 'tot_hd/F')
    
    

    ## Loop over all files
    do_nothing = 0
    last_time = time.time()
    
    for f in file_list:
        data_file = open('%s/%s/%s' % (data_path, cond, f))
        lines = data_file.readlines()
        tmpcnt = 0
        
        # loop over each line
        for line in lines: 
        
            # ignore blank lines
            if len(line.strip()) == 0:
                continue
            
            # ignore file headers
            if tmpcnt < headercnt:
                tmpcnt += 1
                continue
            
            # process frame
            if "#" in line :
                hist3.Fill(linecnt)
                linecnt = 0
                nFrames += 1
        
                # find clusters
                singles, all_clusters = MakeClusters(X, Y, TOT, TOT_HD, hist0, hist1, hist2)
                everyPixel = AllPixels(X, Y, TOT, TOT_HD)
        
                # fill trees
                for cluster in singles: 
                    size1cnt += 1 
                    xt[0] = cluster[0]
                    yt[0] = cluster[1]
                    tott[0] = cluster[2]
                    tott_hd[0] = cluster[3]
                    csize[0] = cluster[4]
                    singleHits.Fill() 
        
                for cluster in all_clusters: 
                    allcnt += 1
                    xt_all[0] = cluster[0]
                    yt_all[0] = cluster[1]
                    tott_all[0] = cluster[2]
                    tott_hd_all[0] = cluster[3]
                    csize_all[0] = cluster[4]
                    allHits.Fill()
        
                for Pixel in everyPixel:
                    everycnt += 1
                    xt_every[0] = Pixel[0]
                    yt_every[0] = Pixel[1]
                    tott_every[0] = Pixel[2]
                    tott_hd_every[0] = Pixel[3]
                    everyPix.Fill()
        
                nclusters = Counting_Clusters(X, Y, TOT, TOT_HD)
                pixel_1hit += nclusters.count(1)
                pixel_2hit += nclusters.count(2)
                pixel_3hit += nclusters.count(3)
                pixel_4hit += nclusters.count(4)
        
                # print progress
                if nFrames%1000 == 0:
                    print "Processed Frame %i (%.5fs/frame)" % (nFrames, (time.time()-last_time)/1000.)
                    last_time = time.time()
                
                # reset    
                X = []
                Y = []
                TOT = []
                TOT_HD = []
        
            # get data from frame
            else : 
                data = line.split()
                X.append(int(data[0]))
                Y.append(int(data[1]))
                TOT.append(float(data[2]))
                TOT_HD.append(float(data[3]))
                linecnt += 1
                del data
           
        del lines 
        data_file.close()

    # print histogrammes
    c = TCanvas("c", "", 0, 0, 800, 700)
    gPad.SetGrid()
    hist0.Draw('')
    hist0.SetLineColor(kRed)
    hist0.SetLineWidth(2) 
    hist0.GetXaxis().SetTitle("Clustersize")
    hist0.GetYaxis().SetTitle("Events")
    hist0.GetYaxis().SetTitleOffset(1.6)
    c.SaveAs('%s/clustersize.pdf' % data_path)

    c1 = TCanvas()
    hist3.Draw()
    hist3.SetLineColor(kRed)
    hist3.SetLineWidth(2) 
    hist3.GetXaxis().SetTitle("Clustersize")
    hist3.GetYaxis().SetTitle("Events")
    hist3.GetYaxis().SetTitleOffset(1.6)
    c1.SaveAs('%s/hitspframe.pdf' % data_path)
        
    c2 = TCanvas("c2", "", 0, 0, 800, 700)
    singleHits.Draw("tot", "tot<40")
    c2.SaveAs('%s/tot_spectrum.pdf' % data_path)
    
    # write to file
    hist0.Write()
    hist1.Write()
    hist2.Write()
    hist3.Write()
    singleHits.Write() 
    allHits.Write()
    everyPix.Write()
    outfile.Close()

    print "fount a total of", everycnt,"pixels and", allcnt, "clusters"
    print "found %i single pixel clusters" % size1cnt
    print "Clusters:", pixel_1hit, "one hit clusters", pixel_2hit, "two hit clusters", pixel_3hit, "three hit clusters", pixel_4hit, "four hit clusters"        



def AllPixels(col, row, tot, tot_hd):
    pixels = [[col[i], row[i]] for i, x in enumerate(col)]
    everyPixel = []
    
    if len(pixels) > 1:
        for hit, TOT, TOT_HD in zip(pixels, tot, tot_hd):
            everyPixel.append([hit[0], hit[1], TOT, TOT_HD])
    else:
        everyPixel = [[pixels[0][0], pixels[0][1], tot[0], tot_hd[0]]]

    return everyPixel



def MakeClusters(col, row, tot, tot_hd, hist0, hist1, hist2):
    oneHitClusters = []
    allHitClusters = []
    pixels = [[col[i], row[i]] for i, x in enumerate(col)]

    #for pixel in pixels :
        #hist1.Fill(pixel[0],pixel[1])

    if len(pixels) > 1:
        results = fclusterdata(pixels, np.sqrt(2.), criterion="distance", method="single") 
        y = np.bincount(results) # histogramm of cluster sizes
        ii = np.nonzero(y)[0]

        do_nothing = 0
        j = 0
        previous = 0

        for result, hit, TOT, TOT_HD in zip(results, pixels, tot, tot_hd):             
            tot_c = 0
            tot_hd_c = 0
            i = 0
            
            # process multi hit clusters
            if y[result] > 1:
                if previous != result:
                    while i <= y[result] - 1:
                        if j < len(results):
                            tot_c += tot[j]
                            tot_hd_c += tot_hd[j]
                            j += 1
                            i += 1
                        if j == len(results):
                             break
                    if tot_c != 0:
                        allHitClusters.append([hit[0], hit[1], tot_c, tot_hd_c, y[result]]) 
                        hist0.Fill(y[result])
            previous = result
            
            # process single hit clusters
            if y[result] == 1:
                if j < len(results):
                    hist0.Fill(y[result])
                    oneHitClusters.append([hit[0], hit[1], TOT, TOT_HD, y[result]])
                    allHitClusters.append([hit[0], hit[1], TOT, TOT_HD, y[result]])
                    j += 1

    if len(pixels) == 1:
        oneHitClusters = [[pixels[0][0], pixels[0][1], tot[0], tot_hd[0], 1]]       

    return oneHitClusters, allHitClusters



def Counting_Clusters(col, row, tot, tot_hd):
    sizes = []
    pixels = [[col[i], row[i]] for i, x in enumerate(col)]

    if len(pixels) > 1:
        results = fclusterdata(pixels, np.sqrt(2.), criterion="distance", method="single")                
        y = np.bincount(results)
        ii = np.nonzero(y)[0]

        j = 0
        previous = 0

        for result, hit, TOT, TOT_HD in zip(results, pixels, tot, tot_hd):
            i = 0
            
            if y[result] > 1:
                if previous != result :  
                    while i <= y[result]-1:
                        if j < len(results):
                            j+=1
                            i+=1
                        if j == len(results):
                             break
                    else :
                         sizes.append(y[result])
            previous = result       
            
            if y[result] == 1:
                if j < len(results) :
                    sizes.append(y[result])
                    j += 1 
    
    else :
        oneHitClusters = [[pixels[0][0], pixels[0][1], tot[0],  tot_hd[0]]]

    return sizes


# Run
ReadFile(assembly, source)


