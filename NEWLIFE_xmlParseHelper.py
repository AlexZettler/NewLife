import re
#import imaplib


def timeParserGetSeconds(inp):

    sec = re.search(r"(\d+)s", inp)
    min = re.search(r"(\d+)m", inp)
    hour = re.search(r"(\d+)h", inp)
    day = re.search(r"(\d+)d", inp)

    counter=0
    #Second
    try:
        counter += int(sec.group(1))
    except:
        pass

    #Minute
    try:
        counter += int(min.group(1))*60
    except:
        pass

    #Hour
    try:
        counter += int(hour.group(1))*60*60
    except:
        pass

    #Day
    try:
        counter += int(day.group(1))*60*60*24
        print("day")
    except:
        pass

    return counter


print(timeParserGetSeconds("1d 2h 1m 30s"))