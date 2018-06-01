import os, sys, glob, json

def get_fls(pth, suf):
    # cur = os.getcwd()
    # os.chdir(pth)
    # abspath=os.path.abspath(pth)
    fls = []
    # for fl in glob.glob(suf, recursive=True): //not help 
        # fls.append(abspath+"/"+fl)
    for root, dirs, files in os.walk(pth):
        for fn in files:
            if fn.endswith(suf):
                 fls.append(os.path.join(root, fn))
    # os.chdir(cur)
    return fls

def write(k, v, o_fl):

    if os.path.isfile(o_fl):
        f = open(o_fl, "r+")
        data = json.load(f)
        f.seek(0, 0)
    else:
        f = open(o_fl, "w")
        data = {}
        
    # if len(v) > 1:
        # v_l = [e for e in v]
    # elif len(v) == 1:
        # v_l = v[0]
    # else:
        # return
    if k not in data:
        data[k] =  []
        # if type(data[k]) == list:
    # else:
        # ov = data[k] # work when data[k] is empty?
        # for e in v:
            # data[k].append(e)
    for e in v:
        data[k].append(e)
    f.truncate()
    json.dump(data, f, indent = 2)
    f.write('\n')
    f.close()



if __name__ == "__main__":
    if len(sys.argv) < 4:
        print ("mkcfg wd DIR SUFFIX KEY")
        print ("mkcfg add KEY MAP")
    else:
        d = sys.argv[1]
        suf = "."+sys.argv[2]
        k = sys.argv[3]
        o = "config.json"
        write(k, get_fls(d, suf), o)

