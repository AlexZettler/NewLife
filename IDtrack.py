#A dictionary object that allows for a value to be stored and returns it's key.

class id_tracker(dict):

    def __init__(self,greatestKey=0):
        super().__init__()

        self.unusedKeys = []
        self.greatestKey = greatestKey

    #stores a value at a key and returns the key given
    def store_val(self, val):
        if isinstance(val, tracked_object):
            if len(self.unusedKeys) != 0:
                key = self.unusedKeys.pop(0)

            elif len(self.unusedKeys) == 0:
                key = self.greatestKey
                self.greatestKey += 1
            else:
                raise KeyError

            tracked_object.__init__(val, key)

            self.__setitem__(key,val)

            return key
        else:
            print("Value was not stored, it was of type {}".format(type(val)))
            raise ValueError

    #removes the value at the given key and configures this key to be reused
    def del_val(self, key):

        if key in self.keys():
            del(self[key])
            self.unusedKeys.append(key)
        else:
            raise ValueError


class tracked_object(object):
    def __init__(self, id: int):
        self.__id = id

    def get_id(self):
        return self.__id

    def validate_id(self):
        return type(self.__id) is int
