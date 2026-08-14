import click
import os
from tqdm import tqdm

@click.command()
@click.option('--path',
                type=click.Path(exists=True),
                help="Path from which find symlinks",
                required=True)

def run (path):
    folders1 = os.listdir(path)

    for folder1 in tqdm(folders1):
        new_path1 = path+folder1
        if os.path.islink(new_path1):
            real_path = os.path.realpath(new_path1)
            os.remove(new_path1)
            os.rename(real_path,new_path1)
            #print(real_path,new_path1)
        else:
            folders2 = os.listdir(new_path1)
            for folder2 in folders2:
                new_path2 = new_path1 +'/'+folder2
                if os.path.islink(new_path2):
                    real_path = os.path.realpath(new_path2)
                    os.remove(new_path2)
                    os.rename(real_path,new_path2)
                    #print(real_path,new_path2)
                else:
                    folders3 = os.listdir(new_path2)
                    for folder3 in folders3:
                        new_path3 = new_path2 +'/'+folder3
                        if os.path.islink(new_path3):
                            real_path = os.path.realpath(new_path3)
                            os.remove(new_path3)
                            os.rename(real_path,new_path3)
                            #print(real_path,new_path3)
                        else:
                            if os.path.isdir(new_path3):
                                folders4 = os.listdir(new_path3)
                                for folder4 in folders4:
                                    new_path4 = new_path3 +'/'+folder4
                                    if os.path.islink(new_path4):
                                        real_path = os.path.realpath(new_path4)
                                        os.remove(new_path4)
                                        os.rename(real_path,new_path4)
                                        #print(real_path,new_path4)                                               

if __name__ == '__main__':
    run()
