import shutil


def cleanx(complexName):
    shutil.rmtree(f'{complexName}/02_frames')
