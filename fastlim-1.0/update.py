
__author__ = "Michele Papucci <mpapucci@umich.edu>"
__version__ = "0.3"


import sys, os, logging, signal, tarfile, time

if sys.version_info >= (3,):
    import urllib.request as urllib2
    import urllib.parse as urlparse
else:
    import urllib2
    import urlparse

def download_file(url, outpath):
    """ download a file from a given url. """

    u = urllib2.urlopen(url)

    scheme, netloc, path, query, fragment = urlparse.urlsplit(url)
    filename = os.path.basename(path)
    if not filename:
        filename = 'temp.tar.gz'

    with open(os.path.join(outpath,filename), 'wb') as f:
        meta = u.info()
        meta_func = meta.getheaders if hasattr(meta, 'getheaders') else meta.get_all
        meta_length = meta_func("Content-Length")
        file_size = None
        if meta_length:
            file_size = int(meta_length[0])
        print("Downloading: {0} Bytes: {1}".format(url, file_size))

        file_size_dl = 0
        block_sz = 8192
        while True:
            buffer = u.read(block_sz)
            if not buffer:
                break

            file_size_dl += len(buffer)
            f.write(buffer)

            status = "{0:16}".format(file_size_dl)
            if file_size:
                status += "   [{0:6.2f}%]".format(file_size_dl * 100 / file_size)
            status += chr(13)
            print status,
        print("Done.                                                           ")

    return filename


def check_install_update(options):
    """ check if fastlim is the most recent version 
    and download the updates when available. """

    # init
    basepath = options['fastlimdir']
    baseurl = 'http://cern.ch/fastlim/'
    oneday = 24 * 3600
    update_interval =  int(options['update_interval']) * oneday
    time_out = 5
    info = {}
    if os.path.exists(os.path.join(basepath,'config','.autoupdate')): #based on MG5 trick
        for tmp in open(os.path.join(basepath,'config','.autoupdate')):
            line = tmp.split()
            if((len(line) >2) and (line[1] == '=')):
                info[line[0]] = int(line[2])
    if 'code_version_major' not in info:
        info['code_version_major'] = options['version_major']
    if 'code_version_minor' not in info:
        info['code_version_minor'] = options['version_minor']
    if 'previous_check' not in info:
        info['previous_check'] = int(time.time()) - update_interval - oneday

    if (int(time.time()) - info['previous_check']) < update_interval:
        return


    logging.info('Checking if Fastlim is up-to-date...')
    
    to_update = 0

    try:
        response = urllib2.urlopen(baseurl + 'fastlim_vers', timeout = time_out)
        website_versions = {}
        lines = response.read().split('\n')
        [ website_versions.update({item[0]: item[2]}) for item in \
            [ elem.split() for elem in lines ] if ((len(item)>2) and (item[1] == '=')) ]  
    except urllib2.URLError:
        logging.error('failed to check server versions')
        # wait 24h before next check
        fout = open(os.path.join(basepath,'config','.autoupdate'),'w')
        fout.write("code_version_major  = %s\n" % info['code_version_major'])
        fout.write("code_version_minor  = %s\n" % info['code_version_minor'])
        fout.write("previous_check   = %s\n" % (int(time.time()) \
            - update_interval + oneday ))
        fout.close()
        logging.info('Will check again in 1 days.')
        return
    

    updated = False
    if ((int(website_versions['code_version_major']) > int(info['code_version_major'])) or \
        (int(website_versions['code_version_minor']) > int(info['code_version_minor']))):
        logging.info('A new version of Fastlim is available!')
        valid = {"yes":True,   "y":True,  "ye":True, "no":False,     "n":False}
        print "Fastlim version %(code_version_major)s.%(code_version_minor)s is available" % website_versions
        while True:
            # print "Do you want to install the update now? [y/n, default=yes]: "
            print "Do you want to install the update now? [y/n]: "
            choice = raw_input().lower()
            # if choice == '' :
            #     choice = 'yes'
            if (choice in valid):
                break
            else:
                # print "Please respond with 'yes' or 'no' (or 'y' or 'n') or press enter."
                print "Please respond with 'yes' or 'no' (or 'y' or 'n')."
        
        if (choice[0] == 'y'):        
            logging.info('Starting update process.')
            if not os.path.exists(os.path.join(basepath,'downloads')):
                os.mkdir(os.path.join(basepath,'downloads'))            
            name = download_file(baseurl + \
                "downloads/fastlim-%(code_version_major)s.%(code_version_minor)s.tar.gz" % website_versions, \
                os.path.join(basepath,'downloads'))
            os.chdir(os.path.join(basepath, '..'))
            tar = tarfile.open(os.path.join(basepath,'downloads',name))
            tar.extractall()
            tar.close()
            os.remove(os.path.join(basepath,'downloads',name))
            updated = True
            info['code_version_major'] = website_versions['code_version_major']
            info['code_version_minor'] = website_versions['code_version_minor']
    else:
        logging.info('No new versions of Fastlim are available')
    

    # update .autopudate with new nubmers
    fout = open(os.path.join(basepath,'config','.autoupdate'),'w')
    fout.write("code_version_major  = %s\n" % info['code_version_major'])
    fout.write("code_version_minor  = %s\n" % info['code_version_minor'])
    fout.write("previous_check  = %s\n" % int(time.time()))
    fout.close()

    # we need to implement a consistency check
    if updated:
        logging.info('A new version of the code has been installed, please relaunch Fastlim')
        sys.exit(0)

    return


