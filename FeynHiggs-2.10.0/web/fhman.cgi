#! /bin/sh
# a trivial cgi script to display man pages
# this file is part of FeynHiggs
# last modified 9 May 05 th

debug=0

if [ $debug = 1 ] ; then
  HOME=/home/pcl301/hahn/fh-22
  PATH=/home/pcl301/hahn/bin/i586-linux:$HOME/build:$PATH
  FH=http://wwwth.mppmu.mpg.de/members/hahn/cgi-bin/fhucc.cgi
else
  HOME=/home/pcl305/members/heinemey/feynhiggs
  PATH=$HOME/bin:$PATH
fi

eval `unescape $QUERY_STRING`
MANFILE=$HOME/man/man1/$man.1

cat - << _EOF_
Content-type: text/html

<html>
<head>
<title>$man</title>
</head>
<body bgcolor=#ffffff>
<h1>$man</h1>
<pre>
_EOF_

/usr/bin/man -l $MANFILE 2>/dev/null | sed "
	/^[A-Z]/d
	/^\o033\[1m/{
	  s:\o033\[1m:</pre><h2>:
	  s:\o033\[0m:</h2><pre>:
	}
	s:\o033\[1m:<b>:g
	s:\o033\[0m:</b>:g
	s:\o033\[4m:<i>:g
	s:\o033\[24m:</i>:g"

cat << _EOF_
</pre>
</body>
</html>
_EOF_

