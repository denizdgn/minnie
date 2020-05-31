import os
import pandas as pd
import logging
import sys
import interfacea as ia
import itertools
import glob
import pathos

# Setup logger
# _private name to prevent collision/confusion with parent logger
logging.getLogger(__name__).addHandler(logging.NullHandler())


def split_pdbs(target,complexName):

    logging.info(f'Splitting {complexName.upper()} ..')
    pathx=os.getcwd()
    if not os.path.exists(f'{complexName}/02_frames'):
        os.makedirs(f'{complexName}/02_frames', exist_ok=True)
    pathxx =  f'{pathx}/{complexName}/02_frames'


    fp = open(target, 'r')
    xx = []
    i = 0
    for line in fp:
        if line.startswith("ATOM"):
            xx.append(line.rstrip())
        elif line.startswith("TER"):
            xx.append(line.rstrip())
        elif line.startswith("HETATM"):
            xx.append(line.rstrip())
        elif line.startswith("END"):
            xx.append("END")
            aa = pd.DataFrame(xx)
            aa.to_csv(pathxx + "/md_" + str(i) + ".pdb", index=None, header=False)
            xx = []
            i = i + 1
            logging.info(f'END stated {i} number of times')

    fp.close()


def paste(x: 'donor_resnm', y: 'donor_resid', a=None, b=None,sep=""):
    """

    :param x: donor_resnm
    :param y: donor_resid
    :param a: acceptor_resnm
    :param b: acceptor_resid
    :return: pasted variables as string
    Ex: paste("a","b") // out: 'ab'
    """
    if a != None :
        return str (x) + str (y) + str(sep) + str (a) + str (b)
    else :
        return str (x)+ str(sep) + str (y)



def comb_int(pdb: "pdb file",complexName,intType:"interaction type",includeIntra,hbond_distance):

    pathx=os.getcwd()

    if not os.path.exists(f'{complexName}/03_interfacea_results/{intType}'):
        os.makedirs(f'{complexName}/03_interfacea_results/{intType}', exist_ok=True)
    pwrite =  f'{pathx}/{complexName}/03_interfacea_results/{intType}'

    ia.set_log_level ('minimal')
    mol = ia.read (pdb)
    analyzer = ia.InteractionAnalyzer (mol)

    if not str(intType) == "hbonds":
        y="get_"+str(intType)
        b= getattr(analyzer,y)
        if str(includeIntra) == "False":
            b()
        else:
            b(include_intra=includeIntra)
        bb = analyzer.itable._table
    else:
        if str(includeIntra) == "False":
            analyzer.get_hbonds (max_distance=float (hbond_distance))
        else:
            analyzer.get_hbonds(include_intra=includeIntra, max_distance=float(hbond_distance))
        bb = analyzer.itable._table

    bonds = bb
    bonds.columns = ['itype', 'donor_chain', 'acceptor_chain', 'donor_resnm', 'acceptor_resnm', 'donor_resid',
                     'acceptor_resid', 'donor_atom', 'acceptor_atom']
    df_table = bonds

    try:
        donor_list=df_table.apply(lambda x: x['donor_resnm'] + str(x['donor_resid']), axis=1)
        df_table.loc[:, 'donor'] = donor_list


        donorC_list=df_table.apply(lambda x,sep="_": x['donor'] +"_"+ str(x['donor_chain']), axis=1)
        df_table.loc[:, 'donorC'] = donorC_list

        acceptor_list=df_table.apply(lambda x: x['acceptor_resnm'] + str(x['acceptor_resid']), axis=1)
        df_table.loc[:, 'acceptor'] = acceptor_list

        acceptorC_list=df_table.apply(lambda x,sep="_": x['acceptor'] +"_"+ str(x['acceptor_chain']), axis=1)
        df_table.loc[:, 'acceptorC'] = acceptorC_list

        donor_acceptor_list = df_table.apply(lambda x, sep=":": x['donorC'] + sep + str(x['acceptorC']), axis=1)
        df_table.loc[:, 'donor_acceptor'] = donor_acceptor_list



        chain_type = df_table.apply(lambda x: "intra" if (x["acceptor_chain"] == x["donor_chain"]) else "inter", axis=1)
        df_table.loc[:, 'chain_type'] = chain_type



        proteinResidues = ['ALA', 'ASN', 'CYS', 'GLU', 'HIS', 'LEU', 'MET', 'PRO', 'THR', 'TYR', 'ARG', 'ASP', 'GLN',
                           'GLY', 'ILE', 'LYS', 'PHE', 'SER', 'TRP', 'VAL']
        dnaResidues = ['DA', 'DC', 'DG', 'DT']

        prot_or_dna = []
        for i in range(len(df_table.donor_resnm)):
            if (((df_table.donor_resnm[i] in dnaResidues) and (df_table.acceptor_resnm[i] in proteinResidues)) or
                    ((df_table.donor_resnm[i] in proteinResidues) and (df_table.acceptor_resnm[i] in dnaResidues))):
                prot_or_dna.append("protein-dna")
            elif (((df_table.donor_resnm[i] in proteinResidues) and (df_table.acceptor_resnm[i] in proteinResidues)) or
                  ((df_table.donor_resnm[i] in proteinResidues) and (df_table.acceptor_resnm[i] in proteinResidues))):
                prot_or_dna.append("protein-protein")
            elif (((df_table.donor_resnm[i] in dnaResidues) and (df_table.acceptor_resnm[i] in dnaResidues)) or
                  ((df_table.donor_resnm[i] in dnaResidues) and (df_table.acceptor_resnm[i] in dnaResidues))):
                prot_or_dna.append("dna-dna")
        df_table["prot_or_dna"] = prot_or_dna




        # renumber index of hbond dataframe
        df_table.index = range (len (df_table))

        if str(intType) == "hbonds" :
            non_spp_atoms = ["O2P", "O1P", "N", "O", "OC1", "OC2", "O4'", "O5'", "O3'", "H", "HA"]
            x = []
            for i in range (len (df_table)) :
                if (df_table["acceptor_atom"][i] or df_table["donor_atom"][i]) in non_spp_atoms :
                    x.append ("non-specific")
                else :
                    x.append ("specific")
            df_table.loc[:, 'specificity'] = x

        time = pdb.split("md_")[1].split(".")[0]
        times = list(itertools.repeat(time, len(df_table)))
        df_table.loc[:, "time"] = times
        logging.info("Writing %s bonds to files...", intType)
        df_table.to_csv(f'{pwrite}/md_{time}.pdb_{intType}_all.csv', index=False)

    except ValueError:
        logging.info(f'Found 0 {intType} interactions !!')


def combine_interfacea_results(complexName):
    pathx=os.getcwd()


    folders = glob.glob(f'{pathx}/{complexName}/03_interfacea_results/*/')
    for dolf in folders:
        bondtype = dolf.split("03_interfacea_results/")[1].split("/")[0]
        files = glob.glob(dolf + "/*all.csv")
        dfx = pd.DataFrame()
        logging.info(f'processing files in {dolf} ..')
        for file in files:
            df = pd.read_csv(file)
            dfx = dfx.append(df)
        logging.info(f'writing to folder {dolf} ..')
        dfx.to_csv(f'{dolf}/{complexName}_merged_{bondtype}.csv', index=False)
        os.chdir(pathx)


    del_files = glob.glob(f'{pathx}/{complexName}/03_interfacea_results/*/*all.csv')
    pool = pathos.multiprocessing.ProcessingPool(pathos.multiprocessing.cpu_count() - 2)
    pool.map(os.remove, del_files)
    pool.close()
