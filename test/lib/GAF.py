import sys

def view_bar(now,total):
    rate = float(now) / total
    rate_win = int(rate*100)
    draw = "\r[%s%s]%d%%" % ("="*int(rate*20)," "*(20-int(rate*20)),rate_win)
    sys.stdout.write(draw)
    sys.stdout.flush()


