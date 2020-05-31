# Copyright 2020 Deniz Dogan, Ezgi Karaca
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Filtering functions.
"""

import pandas as pd
import sys
import glob
import os
import re
import numpy as np
import logging

# Setup logger
# _private name to prevent collision/confusion with parent logger
logging.getLogger(__name__).addHandler(logging.NullHandler())


#inside pathx (MD)
def time_freq_filter(filex,complexName,per):

    pathx = os.getcwd()
    file = os.path.basename(filex)
    fName = complexName
    bondtype = file.split(".csv")[0].split("_merged_")[1]
    first = pd.read_csv(filex)



    os.chdir(pathx)
    if not os.path.exists(f'{complexName}/04_time_freq_filter'):
        os.makedirs(f'{complexName}/04_time_freq_filter', exist_ok=True)
    pathxx=f'{pathx}/{complexName}/04_time_freq_filter'
    os.chdir(pathxx)

    pathy=pathxx+"/"+str(per)+"_freq_filtered"
    if not os.path.exists(str(per)+"_freq_filtered"):
        os.makedirs(str(per)+"_freq_filtered", exist_ok=True)
    os.chdir(pathy)




    if first.empty:
        pathz = pathy + "/" + str(per) + "_freq"
        if not os.path.exists(str(per) + "_freq"):
            os.makedirs(str(per) + "_freq")
        os.chdir(pathz)

        morefirstxy = pd.DataFrame(columns=["donor_acceptor","NumSpp","total","percentage"])
        morefirstxy.to_csv (pathz+"/"+fName+"_"+bondtype+"_"+str(per)+"_freq.csv", index=None)
        os.chdir("..")

        if not os.path.exists(str(per)+"_freq_perres"):
            os.makedirs(str(per)+"_freq_perres")
        pathq=pathy+"/"+str(per)+"_freq_perres"
        os.chdir(pathq)

        first_perres=pd.DataFrame(columns=['itype', 'donor_chain', 'acceptor_chain', 'donor_resnm', 'acceptor_resnm',
                                           'donor_resid','acceptor_resid', 'donor_atom', 'acceptor_atom','chain_type',
                                           "prot_or_dna",'specificity',"time"])
        first_perres.to_csv (pathq+"/"+fName+"_"+bondtype+"_"+str(per)+"_freq_perres.csv", index=None)
    else:

        #fIRST
        logging.info('Finding percentages: {}'.format(fName))
        firstx = []
        for adx in first.donor_acceptor.unique () :
            bbx = first[first["donor_acceptor"] == adx]
            firstx.append([adx,
                            bbx.time.unique().size/first.time.unique().size*100])

        firstxy = pd.DataFrame(firstx)
        firstxy.columns = ["donor_acceptor","percentage"]


        logging.info('Writing to file percentage: {}'.format(fName))
        morefirstxy = firstxy[firstxy.percentage > float(per)]


        if len(morefirstxy.donor_acceptor) == 0:
            pathz = pathy + "/" + str(per) + "_freq"
            if not os.path.exists(str(per) + "_freq"):
                os.makedirs(str(per) + "_freq")
            os.chdir(pathz)

            morefirstxy = pd.DataFrame(columns=firstxy.columns)
            morefirstxy.to_csv (pathz+"/"+fName+"_"+bondtype+"_"+str(per)+"_freq.csv", index=None)
            os.chdir("..")

            if not os.path.exists(str(per) + "_freq_perres"):
                os.makedirs(str(per) + "_freq_perres")
            pathq = pathy + "/" + str(per) + "_freq_perres"
            os.chdir(pathq)

            first_perres= pd.DataFrame(columns=first.columns)
            first_perres.to_csv(pathq + "/" + fName + "_" + bondtype + "_" + str(per) + "_freq_perres.csv", index=None)
        else:

            pathz = pathy + "/" + str(per) + "_freq"
            if not os.path.exists(str(per) + "_freq"):
                os.makedirs(str(per) + "_freq")
            os.chdir(pathz)

            morefirstxy.to_csv (pathz+"/"+fName+"_"+bondtype+"_"+str(per)+"_freq.csv", index=None)
            logging.info('Writing to file list: {}'.format(fName))

            first_perres = pd.DataFrame()
            for da in  morefirstxy.donor_acceptor.unique():
                df = first[first.donor_acceptor == da]
                first_perres=first_perres.append(df)

            first_perres.sort_values(by="time",inplace=True)
            first_perres.reset_index(drop=True)

            os.chdir("..")

            if not os.path.exists(str(per)+"_freq_perres"):
                os.makedirs(str(per)+"_freq_perres")
            pathq=pathy+"/"+str(per)+"_freq_perres"
            os.chdir(pathq)

            first_perres.to_csv (pathq+"/"+fName+"_"+bondtype+"_"+str(per)+"_freq_perres.csv", index=None)


def make_freq_folders(pathy,per):
    """
    Creates folders to write and read common and complex-specific bonds within 05_compare_cx_spp folder
    :param pathy: path to 05_compare_cx_spp
    :param per: time percentage
    """
    import os
    os.chdir(pathy)
    pathz=pathy+"/"+str(per)+"_freq_filtered"
    if not os.path.exists(str(per)+"_freq_filtered"):
        os.makedirs(str(per)+"_freq_filtered",exist_ok=True)



    for fold in ["_freq","_freq_perres"]:
        os.chdir(pathz)
        #to add freq
        pathq=pathz+"/"+str(per)+fold
        if not os.path.exists(str(per)+fold):
            os.makedirs(str(per)+fold,exist_ok=True)
        os.chdir(pathq)

        pathq_common=pathq+"/common"
        if not os.path.exists("common"):
            os.makedirs("common",exist_ok=True)

        os.chdir(pathq)
        pathq_spp=pathq+"/complex_specific"
        if not os.path.exists("complex_specific"):
            os.makedirs("complex_specific",exist_ok=True)


def get_paths(pathy,per,fold,com_spp):
    import os
    os.chdir(pathy)
    PathToWrite = pathy + "/" + per + "_" + "freq_filtered/" + per + fold + "/" + com_spp
    return PathToWrite





def compare_bonds(complexName,per):
    pathx = os.getcwd()
    fName = complexName[0]
    sName = complexName[1]

    file_lists_freq_fName = glob.glob(f'{pathx}/{fName}/04_time_freq_filter/{str(per)}_freq_filtered/{str(per)}_freq/*csv')
    file_lists_freq_sName = glob.glob(f'{pathx}/{sName}/04_time_freq_filter/{str(per)}_freq_filtered/{str(per)}_freq/*csv')

    file_lists_freq = file_lists_freq_fName +  file_lists_freq_sName


    ToCompare = {}

    for filex in file_lists_freq:
        file = os.path.basename(filex)

        if fName in filex:
            Name = fName
        else:
            Name = sName
        bondtype = file.split(f'{Name}_')[1].split("_")[0]
        if bondtype == "ring":
            bondtype = "ring_stacking"
            first = pd.read_csv(filex)
            if bondtype in ToCompare.keys():
                ToCompare[bondtype].update({Name: first})
            else:
                ToCompare.update({bondtype: {Name: first}})

    for bondtype in ToCompare.keys():
        os.chdir(pathx)
        pathy = f'{pathx}/{fName}/05_compare_complex'
        if not os.path.exists(f'{pathx}/{fName}/05_compare_complex'):
            os.makedirs(f'{pathx}/{fName}/05_compare_complex',exist_ok=True)
        os.chdir(pathy)

        pathz = f'{pathx}/{sName}/05_compare_complex'
        if not os.path.exists(f'{pathx}/{sName}/05_compare_complex'):
            os.makedirs(f'{pathx}/{sName}/05_compare_complex',exist_ok=True)
        os.chdir(pathz)

        make_freq_folders(pathy, per)
        fold="_freq"
        morefirstxy = ToCompare[bondtype][fName]

        fold="_freq_perres"
        patha=f'{pathx}/{fName}/04_time_freq_filter/{str(per)}_freq_filtered/{str(per)}{fold}'
        first = pd.read_csv(patha+"/"+fName+"_"+bondtype+"_"+str(per)+fold+".csv")

        #SECOND
        make_freq_folders(pathz, per)
        fold="_freq"
        moresecxy = ToCompare[bondtype][sName]
        logging.info("sName : {}".format(sName))

        fold="_freq_perres"
        patha=f'{pathx}/{sName}/04_time_freq_filter/{str(per)}_freq_filtered/{str(per)}{fold}'
        sec = pd.read_csv(patha+"/"+sName+"_"+bondtype+"_"+str(per)+fold+".csv")

        #find bonds specific to first one
        logging.info("Specific to {}".format(fName))
        i = 0
        spp_first= pd.DataFrame(columns=morefirstxy.columns)
        common_first= pd.DataFrame(columns=morefirstxy.columns)
        for item in morefirstxy.donor_acceptor:
            item_swapped = item.split(":")[1]+":"+item.split(":")[0]
            if item in moresecxy.donor_acceptor.unique():
                common_first = common_first.append(pd.DataFrame(morefirstxy.iloc[i,:]).transpose())
            elif item_swapped in moresecxy.donor_acceptor.unique():
                common_first = common_first.append(pd.DataFrame(morefirstxy.iloc[i,:]).transpose())
            else:
                spp_first = spp_first.append(pd.DataFrame(morefirstxy.iloc[i,:]).transpose())
            i = i+1

        spp_first.sort_values(by="donor_acceptor", ascending=False)
        spp_first.reset_index(drop=True,inplace=True)

        fold="_freq"
        com_spp="complex_specific"
        pathq_spp=get_paths(pathy,str(per),fold,com_spp)
        spp_first.to_csv (pathq_spp+"/"+fName+"_"+bondtype+"_compared_spec.csv", index=False)
        common_first.sort_values(by="donor_acceptor", ascending=False)
        common_first.reset_index(drop=True,inplace=True)

        com_spp="common"
        pathq_common=get_paths(pathy,str(per),fold,com_spp)
        common_first.to_csv (pathq_common+"/"+fName+"_"+bondtype+"_compared_common.csv", index=False)

        #find bonds specific to second one
        logging.info("Specific to  {}".format(sName))
        i = 0
        spp_sec= pd.DataFrame(columns=moresecxy.columns)
        common_sec= pd.DataFrame(columns=moresecxy.columns)
        for item in moresecxy.donor_acceptor:
            item_swapped = item.split(":")[1] + ":" + item.split(":")[0]
            if item in morefirstxy.donor_acceptor.unique():
                common_sec = common_sec.append(pd.DataFrame(moresecxy.iloc[i,:]).transpose())
            elif item_swapped in morefirstxy.donor_acceptor.unique():
                common_sec = common_sec.append(pd.DataFrame(moresecxy.iloc[i,:]).transpose())
            else:
                spp_sec = spp_sec.append(pd.DataFrame(moresecxy.iloc[i,:]).transpose())
            i = i+1
        spp_sec.sort_values(by="donor_acceptor", ascending=False)
        spp_sec.reset_index(drop=True,inplace=True)
        fold="_freq"
        com_spp="complex_specific"
        pathq_spp=get_paths(pathz,str(per),fold,com_spp)
        spp_sec.to_csv (pathq_spp+"/"+sName+"_"+bondtype+"_compared_spec.csv", index=False)
        com_spp = "common"
        pathq_common = get_paths(pathz, str(per), fold, com_spp)
        common_sec.sort_values(by="donor_acceptor", ascending=False)
        common_sec.reset_index(drop=True,inplace=True)
        common_sec.to_csv (pathq_common+"/"+sName+"_"+bondtype+"_compared_common.csv", index=False)


        #find bonds specific to first one
        logging.info("Specific list to {}".format(fName))

        spp_first_perres= pd.DataFrame(columns=first.columns)

        for item in spp_first.donor_acceptor.unique():
            f = first[first.donor_acceptor == item]
            spp_first_perres = spp_first_perres.append(f)


        spp_first_perres.sort_values(by="time", ascending=False)
        spp_first_perres.reset_index(drop=True,inplace=True)

        fold="_freq_perres"
        com_spp="complex_specific"
        pathl_spp=get_paths(pathy,str(per),fold,com_spp)
        spp_first_perres.to_csv (pathl_spp+"/"+fName+"_"+bondtype+"_compared_spec_perres.csv", index=False)
        common_first_perres= pd.DataFrame(columns=first.columns)
        for item in common_first.donor_acceptor.unique():
            f = first[first.donor_acceptor == item]
            common_first_perres = common_first_perres.append(f)
        common_first_perres.sort_values(by="time", ascending=False)
        common_first_perres.reset_index(drop=True,inplace=True)

        com_spp="common"
        pathl_common=get_paths(pathy,str(per),fold,com_spp)
        common_first_perres.to_csv (pathl_common+"/"+fName+"_"+bondtype+"_compared_common_perres.csv", index=False)


        #find bonds specific to second one
        spp_sec_perres= pd.DataFrame(columns=sec.columns)

        for item in spp_sec.donor_acceptor.unique():
            f = sec[sec.donor_acceptor == item]
            spp_sec_perres = spp_sec_perres.append(f)
        spp_sec_perres.sort_values(by="time", ascending=False)
        spp_sec_perres.reset_index(drop=True,inplace=True)
        logging.info(f'Writing to {pathl_spp}')
        logging.info(f'Writing to {pathl_spp}/{sName}_{bondtype}_compared_spec_perres.csv')


        fold="_freq_perres"
        com_spp="complex_specific"
        pathl_spp=get_paths(pathz,str(per),fold,com_spp)
        spp_sec_perres.to_csv (pathl_spp+"/"+sName+"_"+bondtype+"_compared_spec_perres.csv", index=False)
        common_sec_perres= pd.DataFrame(columns=sec.columns)
        for item in common_first.donor_acceptor.unique():
            f = sec[sec.donor_acceptor == item]
            common_sec_perres = common_sec_perres.append(f)
        common_sec_perres.sort_values(by="time", ascending=False)
        common_sec_perres.reset_index(drop=True,inplace=True)

        logging.info(f'Writing to {pathl_common}')
        logging.info(f'Writing to {pathl_common}/{sName}_{bondtype}_compared_common_perres.csv')
        com_spp="common"
        pathl_common=get_paths(pathz,str(per),fold,com_spp)
        common_sec_perres.to_csv (pathl_common+"/"+sName+"_"+bondtype+"_compared_common_perres.csv", index=False)
