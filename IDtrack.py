#A dictionary object that allows for a value to be stored and returns it's key.

class IDtracker(dict):

    def __init__(self,greatestKey=0):
        super().__init__()

        self.unusedKeys = []
        self.greatestKey = greatestKey

    #stores a value at a key and returns the key given
    def StoreVal(self,val):

        if len(self.unusedKeys) != 0:
            key = self.unusedKeys.pop(0)

        elif len(self.unusedKeys) == 0:
            key = self.greatestKey
            self.greatestKey += 1
        else:
            raise KeyError

        val.id = key
        self.__setitem__(key,val)
        return key

    #removes the value at the given key and configures this key to be reused
    def delVal(self,key):
        if key in self.keys():
            del(self[key])
            self.unusedKeys.append(key)
        else:
            raise ValueError

    def sliceMe(self,*args):

        # lets get slicing
        #print("start: {}, stop {}".format(key.start, key.stop))

        if type(args[0])==slice:
            print("start: {},stop: {}".format(args[0].start,args[0].stop))
