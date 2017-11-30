#A dictionary object that allows for a value to be stored and returns it's key.

class IDTracker(dict):
    '''
    A dictionary extension used to automatically assign a free key.
    Values stored must be a TrackedObject subclass, so the value can be acquired from the object itself
    '''
    def __init__(self, greatest_key: int= 0):
        super().__init__()

        #A list of keys that were assigned, but have since been deleted
        self.unusedKeys = []
        #A counter used to track the largest key
        self.greatestKey = greatest_key


    def store_value(self, value) -> int:
        '''
        stores a value at the next key available and returns the given key

        :param value: The TrackedObject that is to be stored
        :return: The key given
        '''
        if isinstance(value, TrackedObject):
            if len(self.unusedKeys) != 0:
                key = self.unusedKeys.pop(0)

            elif len(self.unusedKeys) == 0:
                key = self.greatestKey
                self.greatestKey += 1
            else:
                raise KeyError

            #assign an id to the tracked object
            TrackedObject.__init__(value, key)

            self.__setitem__(key, value)

            return key
        else:
            raise ValueError("Value was not stored, it was of type {}".format(type(value)))

    def delete_value(self, key: int) -> None:
        '''
        Removes the value at the given key and configures this key to be reused

        :param key: the key to be removed
        :return: None
        '''

        if key in self.keys():
            del(self[key])
            self.unusedKeys.append(key)
        else:
            raise ValueError


class TrackedObject(object):
    def __init__(self, _id: int):
        print("TrackedObject")
        self.__id = _id

    def get_id(self):
        return self.__id

    def validate_id(self):
        return type(self.__id) is int
