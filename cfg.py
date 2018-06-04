import os, sys, glob, json

def wd_fls(pth, suf):
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

def update_dict(d, k, vs):
    if len(vs):
        if k in d:
            if type(d[k]) == list:
                for v in vs:
                    if v not in d[k]:
                        d[k].append(v)
            elif type(d[k]) == str:
                ov = d[k]
                if ov != "" and ov not in vs: #even though it could be weird ov can be empty
                    vs.append(ov)
                if len(vs) == 1:
                    d[k] = vs[0]
                else:
                    d[k] = vs
        else:
            if len(vs) == 1:
                d[k] = vs[0]
            else:
                d[k] = vs

def load(fl):
    if os.path.isfile(fl):
        f = open(fl, "r+")
        d = json.load(f)
        f.seek(0, 0)
        f.truncate()
    else:
        f = open(fl, "w")
        d = {}

    return [f, d] 

def dump(d, f):
    f.truncate()
    json.dump(d, f, indent = 2)
    f.write('\n')
    f.close()

# def write(k, v, o_fl):
    # if os.path.isfile(o_fl):
        # f = open(o_fl, "r+")
        # data = json.load(f)
        # f.seek(0, 0)
    # else:
        # f = open(o_fl, "w")
        # data = {}
        
    # if k not in data:
        # data[k] =  []
    # for e in v:
        # data[k].append(e)

def help():
    print ("Usage:    config  wd <dir> <suffix> <key>  # add key-value pairs through wildcard") 
    print ("Usage:    config  ad <key> <value>         # add key-value pairs through values") 


if __name__ == "__main__":
    fn = "config.json"
    if len(sys.argv) < 2:
        help()
    elif sys.argv[1] == "wd":
        if len(sys.argv) < 5:
            help()
        else:
            fldr = sys.argv[2]
            suf = "."+sys.argv[3]
            key = sys.argv[4]
            
            nl = wd_fls(fldr, suf)
            
            [f, dj] = load(fn)
            update_dict(dj, key, nl)
            dump(dj, f)
    elif sys.argv[1] == "ad":
        if len(sys.argv) < 4:
            help()
        else:
            key = sys.argv[2]
            val = sys.argv[3]
            nl = [val]
            
            [f, dj] =load(fn)
            update_dict(dj, key, nl); 
            dump(dj, f)
    else:
        help()

