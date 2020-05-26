#!/usr/bin/env python
import os
import pandas as pd
import logging
import sys
import argparse
import core
from  core import analysis
from  core import filtering
from core import graphs
from  core import clean
import glob
import pathos
from pathos.multiprocessing import ProcessingPool as Pool






###### ---- Preparations ---- ######
parser = argparse.ArgumentParser(prog='minnie')
subparsers = parser.add_subparsers(dest="subcommand", title="Commands", metavar="", help="")


def subcommand(parent=subparsers):

    def decorator(func):

        parser = subparsers.add_parser(func.__name__, help= func.__doc__ , conflict_handler='resolve')
        parserg = parser.add_argument_group("parameters to give")
        parserg.set_defaults(func=func)
        return parserg
    return decorator


def option(*args,**kwargs):
    def decorator(parserg):
        parserg.add_argument(*list(args), **kwargs)
        return parserg
    return decorator





###### ---- Subcommands ---- ######
# ---- split pdbs ---- #
@option('--help','-h',action='store_true')
@option('--complexName','-cn',dest="complexName",nargs="*", help="Project ID of your complex(s)")
@option('--pdbs','-p',nargs="*", help="Give trajectory(s) as *.pdb")
@subcommand()
def splitpdbs(self):
    """Split trajectory into single frames"""
    if self.help:
        print(f'\n\033[1mUsage: minnie splitpdbs \n'
              f'                        -cn, --complexName     <string> <string>     \n '
              f'                                               Project ID of your complex(s)\n\n'

              f'                        -p, --pdbs             [<traj.pdb>] [<traj.pdb>]   \n'
              f'                                               Trajectory ID of your complex(s)\033[0m \n\n\n\n'
          
              f'\n\033[1mUsage example:\033[0m\n\n'
              " minnie splitpdbs -cn sox4 sox18 -p sox4.pdb sox18.pdb \n"
              " minnie splitpdbs -cn sox4  -p sox4.pdb  \n")

    elif not self.pdbs:
        print(f'where is pdb??')
    elif not self.complexName:
        print(f'Please specify complex name(s)')

    elif len(self.pdbs) == 1:
        analysis.split_pdbs(self.pdbs[0], self.complexName[0])
    elif len(self.pdbs) > 1:
        for i in range(len(self.pdbs)):
            analysis.split_pdbs(self.pdbs[i],self.complexName[i])



# ---- find bonds ---- #
@option('--help','-h',action='store_true')
@option('--systematic','-s',choices=["True","False"], help="Use this option if you want to analyze multiple pdbs in a folder",
        default="True")
@option('--pdbs','-p',nargs="*", help="Give single *.pdb or give folder path ")
@option('--complexName','-cn', help="Project ID of your complex(s)")
@option('-i', choices=["hbonds","ionic","hydrophobic","ring_stacking","all"],default="hbonds",
                     dest="intType", help="choose interaction type")
@option('-intra',"--includeIntra", choices=["True","False"],default=False,
                     help="includes intramolecular residue pairs in the search. By default, False.")
@option("-d", default=2.5, dest="hbond_distance",help=" cutoff value to find hbonds")
@subcommand()
def findbonds(self):
    """Calculates interactions between and/or within monomers"""
    if self.help:
        print("Calculates interactions between and/or within monomers\n"
            f'\n\033[1mUsage: minnie findbonds \n'
              f'                        -cn, --complexName     <string>     \n '
              f'                                               Project ID of your complex\n\n'
              
              f'                        -p, --pdbs             [<.pdb>/<path>] (singleframe.pdb)   \n'
              f'                                               Give single *.pdb or give folder path \n\n'
              
              f'                        -i                     [<hbonds>/<ionic>/<hydrophobic>/<ring_stacking>/<all>] (hbonds)    \n'
              f'                                               Calculates which types of interactions \n\n'
              
              f'                        -d                      <float> (2.5)                 \n'
              f'                                               Cut-off to define a hydrogen bond\n\n'
              
              f'                        -intra, --includeIntra [<"True">/<"False">] ("False")  \n'
              f'                                               What do you want to analyze, all or only inter-monomer contacts? \033[0m \n\n\n\n'
                          
              f'\n\033[1mUsage example:\033[0m\n\n'
              " Single frame    - minnie findbonds -cn sox4 -p sox4/02_frames/md_0.pdb -i hbonds  -s False  \n"
              " Multiple frames - minnie findbonds -cn sox4 -p sox4/02_frames/* -i hbonds \n"
              " Multiple frames - minnie findbonds -cn sox4 -p sox4/02_frames/* -i all \n"

              )
    elif not self.pdbs:
        print(f'where is pdb??')
    elif not self.complexName:
        print(f'Please specify complex name(s)')


    elif (self.systematic) == "True":
        pdb_list = self.pdbs
        if (self.intType == "all"):
            for intType in ["hbonds", "ionic", "hydrophobic", "ring_stacking"]:
                pool = Pool(pathos.multiprocessing.cpu_count() - 2)
                pool.map(analysis.comb_int, pdb_list, len(pdb_list) * [str(self.complexName)],
                         len(pdb_list) * [str(intType)], len(pdb_list) * [str(self.includeIntra)],
                         len(pdb_list) * [str(self.hbond_distance)])
                #pool.close()

        else:
            pool = pathos.multiprocessing.ProcessingPool(pathos.multiprocessing.cpu_count() - 2)
            pool.map(analysis.comb_int, pdb_list, len(pdb_list) * [str(self.complexName)],
                     len(pdb_list) * [str(self.intType)],len(pdb_list) * [str(self.includeIntra)],
                     len(pdb_list) * [str(self.hbond_distance)] )
            pool.close()
        analysis.combine_interfacea_results(self.complexName)
    elif (self.systematic) == "False":
        if (self.intType == "all"):
            for intType in ["hbonds", "ionic", "hydrophobic", "ring_stacking"]:
                analysis.comb_int(self.pdbs[0], self.complexName, intType, self.includeIntra, self.hbond_distance)
        else:
            analysis.comb_int(self.pdbs[0],self.complexName,self.intType,self.includeIntra,self.hbond_distance)


        analysis.combine_interfacea_results(self.complexName)





# ---- Time filtering ---- #
@option('--help','-h',action='store_true')
@option('--files','-f',nargs="*", help="Files")
@option('--complexName','-cn',help="Project ID of your complex(s)")
@option('--per', help="observation frequency (in %) to classify an interaction as critical?", type = int)
@subcommand()
def timefilter(self):
    """Apply critical interaction filter"""
    path_back = os.getcwd()
    if self.help:
        print(f'\n\033[1mUsage: minnie timefilter \n'
              f'                        -cn, --complexName     <string>      \n '
              f'                                               Project ID of your complex\n\n'
              
              f'                         -f, --files            [<.csv>]               \n'
              f'                                                Files of your complex\n\n'
              f'                         --per                  <float>                \n'
              f'                                                Observation frequency to classify an interaction as critical\033[0m \n\n\n\n'

            f'\n\033[1mUsage example:\033[0m\n\n'
              " Single file    - minnie timefilter -f sox4/03_interfacea_results/hbonds/sox4_merged_hbonds.csv -cn sox4 --per 25  \n"
              " Multiple files - minnie timefilter -f sox4/03_interfacea_results/*/sox4_merged_*.csv -cn sox4 --per 25  \n"
              )
    elif not self.files:
        print(f'\nwhere is the file(s) ??\n')
    elif not self.complexName:
        print(f'\nPlease specify a complex name(s) !!\n')
    elif not self.per:
        print(f'\nPlease specify a cutoff value to filter out bonds !!\n')

    if (self.per):
        for filex in self.files:
            os.chdir(path_back)
            filtering.time_freq_filter(filex,self.complexName,self.per)
        #pool = pathos.multiprocessing.ProcessingPool(pathos.multiprocessing.cpu_count() - 2)
        #pool.map(filtering.time_freq_filter, self.files,
        #         len(self.files)*[self.complexName],
        #         len(self.files)*[self.per])
        #pool.close()


# ---- Compare interaction networks between two complex ---- #
@option('--help','-h',action='store_true')
@option('--per', help="observation frequency (in %) to classify an interaction as critical?", type = int)
@option('--complexName','-cn',nargs=2, help="Project ID of your complex(s)")
@subcommand()
def compareCX(self):
    """Calculate common and distinct interactions in two cases"""
    if self.help:
        print(f'\n\033[1mUsage: minnie compareCX  \n'
              f'                        -cn, --complexName     <string> <string>     \n '
              f'                                               Project ID of your complex(s)\n\n'
              f'                        --per                  <float>                \n'
              f'                                                Observation frequency to classify an interaction as critical\033[0m \n\n\n\n'
              
              f'\n\033[1mUsage example:\033[0m\n\n'
              " minnie compareCX -cn sox4 sox18 --per 25  \n")
    else:
        filtering.compare_bonds(self.complexName,self.per)





# ---- draw graphs! ---- #
@option('--help','-h',action='store_true')
@option('--per', help="observation frequency (in %) to classify an interaction as critical?", type = int)
@option('--colors',nargs="*", help="Color IDs of the complexes (for plotting)",default=['#D9B4CC', '#6F81A6'])
@option('--chainIDs','-c',nargs=2,dest="chains", help="Give chainID(s)")
@option('--complexName','-cn',nargs=2, help="Project ID of your complex(s)")
@option('-s', choices=["specific","common"],default="specific", dest="spp")
@option('-i', choices=["hbonds","ionic","hydrophobic","ring_stacking","all"],default="hbonds",
                     dest="intType", help="choose interaction type")
@option('-b',"--between", choices=["protein-dna","all"], dest="between",default="all")
@option('--filename', dest="filename",help="filename to use while writing graph",default="")
@subcommand()
def graph(self):
    """aaaaand graphs!"""
    try:
        if self.help:
            print(f'\n\033[1mUsage: minnie graph  \n'
                  f'                        -cn, --complexName     <string> <string>     \n'
                  f'                                               Project IDs of your complex(s)\n\n'
                  
                  f'                        --per                  <float>                \n'
                  f'                                               Observation frequency to classify an interaction as critical \n\n'
    
                  
                  f'                        -b, --between          [<protein-dna>/<all>] (all)   \n'
                  f'                                               Between protein-dna or keep all \n\n'   
                  
                  f'                        -c, --chainIDs         <string> <string>       \n'
                  f'                                               Give ChainIDs to proceed\n\n'
                  
                  f'                        --filename             <string>                           \n'
                  f'                                               Give a name to output file (optional)\n\n'
                  
                  f'                        --colors               [<hex colors>] ("#D9B4CC", "#6F81A6")     \n'
                  f'                                               Color IDs of the complexes (optional)\n\n'
                  
                  
                  f'                        -i                     [<hbonds>/<ionic>/<hydrophobic>/<ring_stacking>/<all>] (hbonds)    \n'
                  f'                                               Calculates which types of interactions \n\n'
    
                  f'                        -s                     [<specific>/<common>] (specific)                           \n'
                  f'                                               Complex-specific or common interactions\033[0m \n\n\n\n'
    
                f'Please do not use "--between" and "--chainIDs" options at the same time\n\n'
    
                f'\n\033[1mUsage example:\033[0m\n\n'
                  " minnie graph -cn 'sox4' 'sox18' --per 25 -i hbonds -s specific  -c A+B C --colors '#D9B4CC' '#6F81A6' \n"
                  " minnie graph -cn 'sox4' 'sox18' --per 25 -i ionic -c A+B C  \n"
                  " minnie graph -cn 'sox4' 'sox18' --per 25 -i ionic -b protein-dna \n"
                  " minnie graph -cn 'sox4' 'sox18' --per 25 -i ionic -b protein-dna --filename sox4_sox18_protein_dna \n")
        elif (self.between) and (self.chains):
            raise Exception()
        elif self.intType == "all":
            for intTypex in ["hbonds","ionic","hydrophobic","ring_stacking"]:
                if  (self.between) == "all":
                    print(intTypex)
                    df_collec=graphs.filter_todraw(self.complexName,self.chains,self.spp,self.per,str(intTypex))
                    graphs.draw_fig(df_collec, str(intTypex), self.complexName[0], self.complexName[1],
                                    self.colors[0], self.colors[1],self.filename, self.spp)

                else:
                    df_collec = graphs.filter_todnaall(self.complexName, self.between, self.spp, self.per, intTypex)
                    graphs.draw_fig(df_collec, intTypex, self.complexName[0], self.complexName[1],
                                    self.colors[0], self.colors[1], self.filename, self.spp)

        elif self.between == "protein-dna":
            df_collec=graphs.filter_todnaall(self.complexName,self.between,self.spp,self.per,self.intType)
            graphs.draw_fig(df_collec, self.intType, self.complexName[0], self.complexName[1],
                            self.colors[0], self.colors[1],self.filename, self.spp)

        else:
            df_collec=graphs.filter_todraw(self.complexName,self.chains,self.spp,self.per,self.intType)
            graphs.draw_fig(df_collec, self.intType, self.complexName[0], self.complexName[1],
                            self.colors[0], self.colors[1],self.filename, self.spp)
    except TypeError:
        print(f'\nPlease check given parameters''')

    except Exception:
        print(f'\nPlease, either use -b or -c option''')









if __name__ == "__main__":
    args = parser.parse_args()
    if args.subcommand is None:
        parser.print_help()
    else:
        args.func(args)



