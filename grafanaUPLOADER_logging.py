#!/usr/bin/env python                                                                                                                                                                                                                                               
# coding=utf-8                                                                                                                                                                                                                                                      
# vim: ts=4 sw=4 et                                                                                                                                                                                                                                                 
# Original Author: Alba Domi - COPYRIGHT 2018                                                                                                                                                                                                                       
# Modified by T. Chiarusi 2018, 2019/04/24 Francescos                                                                                                                                                                                                               
# Compatible with KM3PIPE 9 beta                                                                                                                                                                           

import logging                                                         

logging.basicConfig(filename='/logs/logfile_2.log', filemode='a', format='%(message)s')

import signal
import sys
import subprocess
import re
import collections
from collections import defaultdict
import io
import os
import time
import km3pipe as kp
import numpy as np
from km3pipe.io.daq import TMCHData
from km3pipe import Module
import pandas as pd
import km3db as db
#from km3pipe.core import Pump                                                                                                                                                                                                                                      
from datetime import datetime
#from datetime import datetime as dt                                                                                                                                                                                                                                
from dotenv import load_dotenv


UDP_RATE_PER_DOM = 10  # Hz                                                                                                                                                                                                                                         
LIMIT_TIME_SYNC  = 3    #3701 # .5  # s   value tuned after several trials (3700 trigger the warning printing)                                                                                                                                                     
TIME_OFFSET = 0 # 32.0481 # s                           
write_interval = 30 #time to wait before writing dataframes on csv                                                                                                                                                               
active_doms = []          
                                   

#added_second = datetime.timedelta(0,1)                                                                                                                                                                                                                             
#logging.basicConfig(filename='prova.log',filemode='w',force=True)                                                                                                                  


class UDPAnalyser(Module):

    def configure(self):

      # self.detector = self.require("detector")                                                                                                                                                                                                                    
        self.detector_oid = 4 #db.detectors[db.detectors.SERIALNUMBER == det_id].OID.iloc[0]                                                                                                                                                                        

        self.interval = 30  # seconds of DOWNSAMPLING                                                                                                                                                                                                                            
        self.run_id = defaultdict(int)
        self.timestamp = defaultdict(int)
        self.total_number_udps = defaultdict(int)
        self.data_time = defaultdict(int)
        self.start_of_run_dom = defaultdict(int)
        self.start_of_run_machine = defaultdict(int) #time of the machine when 1st packet arrives --> corresponds to start_of_run_dom, but with the machine time and not with data time                                                                             
        self.end_of_run_dom = defaultdict(int)
        self.end_of_run_machine = defaultdict(int)
        self.run_duration_dom = defaultdict(int)
        self.time = defaultdict(int)
        self.copy_time = defaultdict(int)
        self.time_of_previous_check = defaultdict(int)
        self.first_packet = defaultdict(int)
        self.packet_machine_time = defaultdict(int)
        self.dt_packets = defaultdict(int)
        self.ratio = defaultdict(int)
        self.loss_ratio = defaultdict(int)
        self.timedifference = defaultdict(int)
        self.n_expected_packets = defaultdict(int)
        self.machine_expected = defaultdict(int)
        self.filecounter = 0
        self.dt = defaultdict(int)
        self.datatag = defaultdict(int)
        
        self.timevalid = defaultdict(int) #timevalid bit value and its copy
        self.copy_timevalid = defaultdict(int) 
        

        self.total_missing_udps = defaultdict(int) #total number of missing packets as counted by the 100ms check function
        self.clock_reset_counter = defaultdict(int)

        self.flag_100ms = 0
        self.flag_delay = 0
        self.flag_timevalid = defaultdict(int)

        self.writing_timestamp = 0
               
        #self.filecounter = 0
        #self.duname=int(self.get('duname'))
        self.filename = ""
        self.du_directory=sys.argv[1]
        self.time_run_change = 0
        self.duname=int(re.split('(\d+)',self.du_directory)[1])
        detector=os.getenv('DETECTOR')
        #self.clbmap = db.CLBMap("D_BCI0004")
        self.clbmap = db.CLBMap(detector)
        self.Dom_id_name = 0

        self.testdf=pd.DataFrame(columns=["DU","det_id","run","source","time_valid","tslice_duration","machine_time","packet_time","delay","check100ms","loss_ratio_datatime","loss_ratio_machinetime","total_udp_packets_arrived","expected_from_datatime","expected_from_machinetime","clock_reset_counter"])

    def process(self, blob):

        arrival_time = blob["CHPrefix"].timestamp
        tmch_data = kp.io.daq.TMCHData(io.BytesIO(blob['CHData']))
        dom_id = tmch_data.dom_id
        bit_timevalid = int(bin(tmch_data.dom_status_0)[2]) #bit 31 -> time valid                                                                                                                          
        tslice_duration = 100 #tmch_data.ts_duration_ms #in microseconds  

        #------------------------------------------------------------------------------                                                                                                                

        doms = self.clbmap.dom_ids[dom_id]
        self.Dom_id_name = doms.floor

        if not self.run_id[dom_id]: #fill run_id dictionary with run values;@ 1st timeslice per DOM of the run it enters here (run_id = 0) then id is associated and never enters here until run change                                                             
            self.run_id[dom_id] = tmch_data.run

        if tmch_data.run != self.run_id[dom_id]:
            self.reset_data_end_of_run(dom_id)
            self.run_id[dom_id] = tmch_data.run
            self.start_of_run_dom[dom_id] = (tmch_data.utc_seconds*1e9 + tmch_data.nanoseconds)*1e-6
            self.start_of_run_machine[dom_id] = round(time.time(),1)
            print('NEW RUN',self.run_id[dom_id],' (1st detected packet) for DOM',self.Dom_id_name,' START AT ',self.start_of_run_dom[dom_id], "ms")
            print("MACHINE START TIME = ",self.start_of_run_machine[dom_id])


        if not self.start_of_run_dom[dom_id]:  #for the 1st run, initialize the start time in seconds for each DOM. For subsequent run this is done at run change                                                                                                   
            self.start_of_run_dom[dom_id] = (tmch_data.utc_seconds*1e9 + tmch_data.nanoseconds)*1e-6  #NOTE: 'start of run' actually means the time of 1st analyzed packet. If this program is launched when run is already started, start_of_run is NOT the acual start of run time                                                                                                                                                                                                                                                   
            
            self.start_of_run_machine[dom_id] = round(time.time(),1) #need this for the 1st check interval = [ machine time at start, machine time at start + self.interval ]                                                                                    
             #if self.Dom_id_name==2: #to select a specific dom to print                                                                                                                                                                               

            print('NEW RUN',self.run_id[dom_id],' (1st detected packet) for DOM',self.Dom_id_name,' START AT ',self.start_of_run_dom[dom_id], "ms")
            print("MACHINE START TIME = ",self.start_of_run_machine[dom_id])

        #if self.filename == "":
        #    self.filename = "/home/km3net/analysis/grafanarepo/test_grafana_RUN" + str(tmch_data.run)+ ".csv"
        #    self.write_header(tmch_data.run)
            print(self.filename)

        if not self.timestamp[dom_id]: #enter here at first round and assign value to timestamp[dom_id]; timestamp is increased by 10 s, it's the interval time
            self.reset_data(dom_id)
        #----------------------------------------------------------------                                                                                                                                                                                          
 
        self.total_number_udps[dom_id] += 1 #used to count the total number of udp packets, never reset                                                                                                                                                             
        #----------------------------------------------------------------                                                                                                                                                                                           
        total_time = (tmch_data.utc_seconds*1e9 + tmch_data.nanoseconds)*1e-6  # ms                                                                                                                                                                                 

        self.packet_machine_time[dom_id] = round(time.time(),3)
        self.time[dom_id] = total_time
        
        self.dt[dom_id] = round(self.packet_machine_time[dom_id] - self.time[dom_id]*1e-3,3)

        #----------------------------------------------------------------                                                                                                                                                                                           
        self.timevalid[dom_id] = bit_timevalid
        
        #set flag if time invalid in order to take data from now on until end of run
        if self.timevalid[dom_id] == 0: 
            self.flag_timevalid[dom_id] = 1

        #----------------------------------------------------------------                                                                                                                               
        self.check_packet_loss(dom_id)
        self.check_100ms_sync(arrival_time,dom_id)
        self.check_packets_vs_machine(dom_id)
        self.time_valid_check(dom_id,arrival_time)

        self.copy_time[dom_id] = total_time

        self.copy_timevalid[dom_id] = bit_timevalid

        if self.first_packet[dom_id] == 0:
            self.first_packet[dom_id] = 1
            if self.Dom_id_name not in active_doms:
                active_doms.append(self.Dom_id_name)
                print(active_doms)

        #----------------------------------------------------------------  
        
        if self.writing_timestamp == 0:
            self.writing_timestamp = time.time()

        #----------------------------------------------------------------                                                                                                                                                                                           
        #keep ATTENTION to the difference between data_time(s) and times(ms)                                                                                                                                                                                        
        self.data_time[dom_id] = total_time*1e-3

        #----------------------------------------------------------------                                                                                                                                                                                           
        #function execution > check done for each dom for each TS                                                                                                                                                                                                   
        self.check_data_machine_time(arrival_time,dom_id)
        
        #add data point to dataframe if some kind of unexpected behavior occurred                                                                                                                         
        if self.flag_100ms == 1 or self.flag_delay == 1 or self.flag_timevalid == 1:
            self.testdf.loc[len(self.testdf)+1] = [self.duname,self.detector_oid,self.run_id[dom_id],self.Dom_id_name,bit_timevalid,tslice_duration,self.packet_machine_time[dom_id],self.time[dom_id],self.dt[dom_id],self.dt_packets[dom_id],self.loss_ratio[dom_id],self.ratio[dom_id],self.total_number_udps[dom_id],self.n_expected_packets[dom_id],self.machine_expected[dom_id],self.clock_reset_counter[dom_id]]

        #add data point to dataframe every self.interval, if no unexpected behavior occurred                                                                                                             
        elif self.return_timedelta(dom_id)>self.interval:
            self.testdf.loc[len(self.testdf)+1] = [self.duname,self.detector_oid,self.run_id[dom_id],self.Dom_id_name,bit_timevalid,tslice_duration,self.packet_machine_time[dom_id],self.time[dom_id],self.dt[dom_id],self.dt_packets[dom_id],self.loss_ratio[dom_id],self.ratio[dom_id],self.total_number_udps[dom_id],self.n_expected_packets[dom_id],self.machine_expected[dom_id],self.clock_reset_counter[dom_id]]
            self.reset_data(dom_id)
            
        #write dataframe to csv every write_interval
        if self.return_writing_timedelta() > write_interval:
            self.filename = os.getcwd()+'/grafanarepo/'+self.du_directory+'/test_grafana_RUN'+str(self.run_id[dom_id])+"_"+ str(self.filecounter)+ ".csv"
            print('writing dataframe ',self.filename)
            self.write_data_into_file()
            self.filecounter+=1
            self.testdf = self.testdf.iloc[0:0]
            self.reset_writing_timestamp()

        #reset flags to 0, prepare for next packet                                                                                                                                                       
        self.flag_100ms = 0
        self.flag_delay = 0

        #----------------------------------------------------------------                                                                                                                                                                                           
        return blob

    #----------------------------------------------------------------                                                                                                                                                                                                
    #making the difference of system time.time(s) and timestamp; timestamp changes every self.interval s                                                                                                                                                            
    def return_timedelta(self, dom_id):
        return time.time() - self.timestamp[dom_id]

    #----------------------------------------------------------------                                                                                                                                                                                                        
    def return_writing_timedelta(self):
        return time.time() - self.writing_timestamp
    #----------------------------------------------------------------                                                                                                                                                                                   
    def reset_data(self, dom_id): #called at the end of every self.interval                                                                                                                                                                                         
        self.timestamp[dom_id] = time.time()
        
    #----------------------------------------------------------------                                                                                                                                                                                                        
    def reset_writing_timestamp(self): #called after writing dataframes on csv files
         self.writing_timestamp = 0

    #----------------------------------------------------------------                                                                                                                                                                                               
    def write_header(self, run_id):
        """                                                                                                                                                                                                                                                         
        This function is called once (if the outputfile doesn't exist).                                                                                                                                                                                             
        """
        if not os.path.exists(self.filename):
            out = open(self.filename, "w+")
            out.write("tag\tdet_id\trun\tsource\ttime_valid\ttslice_duration")
            out.write("\tmachine_time\tpacket_time\tdelay\tcheck100ms\ttransmission_ratio_datatime\ttransmission_ratio_machinetime\total_udp_packets_arrived\texpected_from_datatime\texpected_from_machinetime\tclock_reset_counter")
            out.write("\n")
            out.close()

    #----------------------------------------------------------------                                                                                                                                                                                               
    def time_valid_check(self,dom_id,arrival_time):

        """
        Check if a change in the bit timevalid value has occurred                                                                                                                                         
        """
    
        if self.first_packet[dom_id] == 0:
            print("Inside check_timevalid for DOM",self.Dom_id_name,", 1st packet")

        elif self.timevalid[dom_id] != self.copy_timevalid[dom_id]:

                #convert time since epoch to date                                                                                                                                                          
                arrival_string = "{:.3f}".format(arrival_time)
                data_time_string = "{:.3f}".format(round(self.time[dom_id]*1e-3,3))

                arrival_conv = (time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(arrival_time)))+"."+arrival_string.split(".")[1]
                data_time_conv = (time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(round(self.time[dom_id]*1e-3,3))))+"."+data_time_string.split(".")[1]
                logging.error(f'{arrival_conv}, RUN {self.run_id[dom_id]} - DU {self.duname}, DOM {self.Dom_id_name} - BIT TIMEVALID value changed from {self.copy_timevalid[dom_id]} to {self.timevalid[dom_id]}. Packet timestamp = {data_time_conv}')


    #----------------------------------------------------------------                                                                                                                                                                                               
    def check_packet_loss(self, dom_id, end_run=0):
        """                                                                                                                                                                                                                                                         
        Check if the effective number of packets received by each dom                                                                                                                                                                                               
        is different from what expected.                                                                                                                                                                                                                            
        """
        if  self.first_packet[dom_id] == 0:
            print("Inside check_packet_loss for DOM",self.Dom_id_name, ", 1st packet.")

        else:
            self.dt_packets[dom_id] = self.time[dom_id] - self.copy_time[dom_id]
            
            self.n_expected_packets[dom_id] = ((self.time[dom_id]-self.start_of_run_dom[dom_id])/100)+1  #expected number of packets from start time. +1 is for adding 1st starting packet   
            observed_packets = self.total_number_udps[dom_id]
            self.loss_ratio[dom_id] = observed_packets / self.n_expected_packets[dom_id]
            #if self.n_expected_packets[dom_id] != observed_packets:
                #log.warning("PACKET LOSS for DOM {0} = {1} missing packets (this is cumulative)".format(self.Dom_id_name,self.n_expected_packets[dom_id]-observed_packets))
                #print('PACKET LOSS ratio for DOM {0} --> expected {1}, got {2}, ratio is {3}'.format(self.Dom_id_name,self.n_expected_packets[dom_id],observed_packets,self.loss_ratio[dom_id]))              
                #print("PACKET LOSS for DOM {0} = {1} missing packets (this is cumulative)".format(self.Dom_id_name,self.n_expected_packets[dom_id]-observed_packets))
                
    #----------------------------------------------------------------                                                                                                                                                                                               
    def check_100ms_sync(self,arrival_time,dom_id):
        """                                                                                                                                                                                                                                                         
        Check if there are some consecutive udp packets                                                                                                                                                                                                             
        with delta_t != 100 ms.                                                                                                                                                                                                                                     
        """
        if  self.first_packet[dom_id] == 0:
            print("Inside check_100_ms for DOM",self.Dom_id_name,", 1st packet")

        else:
            self.dt_packets[dom_id] = round(self.time[dom_id],3) - round(self.copy_time[dom_id],3)
            if self.dt_packets[dom_id] < -946080000:
                self.clock_reset_counter[dom_id]+=1
                #convert time since epoch to date
                arrival_string = "{:.3f}".format(arrival_time)
                data_time_string = "{:.3f}".format(self.data_time[dom_id])

                arrival_conv = (time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(arrival_time)))+"."+arrival_string.split(".")[1]
                data_time_conv = (time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(self.data_time[dom_id])))+"."+data_time_string.split(".")[1]

                logging.error(f'{arrival_conv} RUN {self.run_id[dom_id]} - DU {self.duname}, DOM {self.Dom_id_name} - START of CLOCK RESET ERROR: packet time = {data_time_conv} -> delay = {self.timedifference[dom_id]}')

            elif self.dt_packets[dom_id] > 946080000:
                
                #convert time since epoch to date
                arrival_string = "{:.3f}".format(arrival_time)
                data_time_string = "{:.3f}".format(self.data_time[dom_id])

                arrival_conv = (time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(arrival_time)))+"."+arrival_string.split(".")[1]
                data_time_conv = (time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(self.data_time[dom_id])))+"."+data_time_string.split(".")[1]

                logging.error(f'{arrival_conv} RUN {self.run_id[dom_id]} - DU {self.duname}, DOM {self.Dom_id_name} - END of CLOCK RESET ERROR: packet time = {data_time_conv} -> delay = {self.timedifference[dom_id]}')


            # print("dt between packets for DOM",self.Dom_id_name," = ",self.dt_packets[dom_id])
            elif self.dt_packets[dom_id] != 100:
                self.flag_100ms = 1 #raise flag to write data point to dataframe                                                                                                                          
                self.total_missing_udps[dom_id]+= (self.dt_packets[dom_id]-100)/100

                #convert time since epoch to date                                                                                                                                                         
                arrival_string = "{:.3f}".format(arrival_time)
                data_time_string = "{:.3f}".format(round(self.time[dom_id]*1e-3,3))
                previous_time_string = "{:.3f}".format(round(self.copy_time[dom_id]*1e-3,3))

                arrival_conv = (time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(arrival_time)))+"."+arrival_string.split(".")[1]
                data_time_conv = (time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(round(self.time[dom_id]*1e-3,3))))+"."+data_time_string.split(".")[1]
                previous_time_conv = (time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(round(self.copy_time[dom_id]*1e-3,3))))+"."+previous_time_string.split(".")[1]

                logging.error(f'{arrival_conv} RUN {self.run_id[dom_id]} - DU {self.duname}, DOM {self.Dom_id_name} - 100ms ERROR: dt = {self.dt_packets[dom_id]} -> current packet timestamp = {data_time_conv}, previous packet timestamp = {previous_time_conv}')


    #----------------------------------------------------------------                                                                                                                                                                                               
    def check_packets_vs_machine(self, dom_id):
        """                                                                                                                                                                 
        Check packets expected from data elapsed time vs
        expected number considering machine elapsed time                                                                                                                                                                                                            
        """
        if  self.first_packet[dom_id] == 0:
            print("Inside check_packet_vs_machine for DOM",self.Dom_id_name,", 1st packet")
        else:
            observed = self.total_number_udps[dom_id]
            self.machine_expected[dom_id] = (round((self.packet_machine_time[dom_id]- self.start_of_run_machine[dom_id]),2)*10)+1
            self.ratio[dom_id] = observed/self.machine_expected[dom_id]
            
    #----------------------------------------------------------------                                                                                                                                                                                               
    def check_data_machine_time(self, arrival_time, dom_id): #arrival time on machine vs time of data                                                                                                                                                               
        """                                                                                                                                                                                                                                                         
        Check if the timestamp of each udp packet and                                                                                                                                                                                                               
        its arrival time on the machine is > 1 minute.                                                                                                                                                                                                              
        """
        self.timedifference[dom_id] = arrival_time - self.data_time[dom_id]+TIME_OFFSET
        #if self.Dom_id_name==1:
        #    print('check clock reset for DOM1, dt = ',self.timedifference[dom_id])
#        if abs(self.timedifference[dom_id]) > 946080000: #if delay is larger than 30 years                                                                                                         
            #convert time since epoch to date                                                                                                                                                             
#            arrival_string = "{:.3f}".format(arrival_time)
#            data_time_string = "{:.3f}".format(self.data_time[dom_id])

#            arrival_conv = (time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(arrival_time)))+"."+arrival_string.split(".")[1]
#            data_time_conv = (time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(self.data_time[dom_id])))+"."+data_time_string.split(".")[1]

#            logging.error(f'{arrival_conv} RUN {self.run_id[dom_id]} - DU {self.duname}, DOM {self.Dom_id_name} - CLOCK RESET ERROR: packet time = {data_time_conv} -> delay = {self.timedifference[dom_id]}')

        if abs(self.timedifference[dom_id]) > LIMIT_TIME_SYNC:
            self.flag_delay = 1

            #convert time since epoch to date                                                                                                                                                            
            arrival_string = "{:.3f}".format(arrival_time)
            data_time_string = "{:.3f}".format(self.data_time[dom_id])

            arrival_conv = (time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(arrival_time)))+"."+arrival_string.split(".")[1]
            data_time_conv = (time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(self.data_time[dom_id])))+"."+data_time_string.split(".")[1]

            logging.error(f'{arrival_conv} RUN {self.run_id[dom_id]} - DU {self.duname}, DOM {self.Dom_id_name} - DELAY ERROR: packet time = {data_time_conv} -> delay = {self.timedifference[dom_id]}')

    #----------------------------------------------------------------                                                                                                                                                                                                   
    def write_data_into_file(self):
        self.testdf.to_csv(self.filename, mode='a',index=False)
        
    #----------------------------------------------------------------                                                                                                                                                                                               
    def reset_data_end_of_run(self, dom_id): #called only at the end of the run                                                                                                                                                                                     
        self.run_id[dom_id] = 0
        self.total_number_udps[dom_id] = 0
        self.total_missing_udps[dom_id] = 0
        self.timestamp[dom_id] = time.time()
        self.first_packet[dom_id] = 0 
        self.start_of_run_dom[dom_id] = 0
        self.start_of_run_machine[dom_id] = 0
        self.filecounter = 0
        self.dt_packets[dom_id] = 0
        self.flag_timevalid[dom_id] = 0
        
    #    self.clock_reset_counter[dom_id] = 0
    #----------------------------------------------------------------                                                                                                                                                                                               
def signal_handler(sig, frame):
        print('You pressed Ctrl+\!')
        print("\nMON quitting time: {0} - {1}\n".format(datetime.timestamp(datetime.now()),datetime.now()));
        sys.exit(0)
    #----------------------------------------------------------------                                                                                                                                                                                               

def main():
    #detector = kp.hardware.Detector(det_id = 29)                                                                                                                                                                                                                   
    load_dotenv()
    host_ip = os.getenv('HOST_IP')
    signal.signal(signal.SIGQUIT, signal_handler)
    print("\n---- Tom's modified udpAnalyser v1.0 ----\n")
    print("\nMON starting time: {0} - {1}".format(datetime.timestamp(datetime.now()),datetime.now()));
    print("\npress Ctrl+\ (SIGQUIT)for quitting gently and getting the stopping time.")
    pipe = kp.Pipeline(timeit=True)
    pipe.attach(kp.io.CHPump, host='jpp',
                              port=9999,
                              tags=sys.argv[1],
                              timeout=60*60*3,
                              max_queue=10000000000)
    pipe.attach(UDPAnalyser)

    pipe.drain()

if __name__ == "__main__":
    main()
