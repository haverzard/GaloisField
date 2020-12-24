from exceptions import FPNegativeDegreeNotAllowed
from util import get_sign


class FastPolynom:
    def __init__(self, container=None):
        self.container = container or {}

    def deepcopy(self):
        return FastPolynom(self.container.copy())

    def broadcast_modulo(self, p):
        for i in self.container:
            self.container[i] %= p

    def get_max_degree(self):
        keys = list(self.container.keys())
        if keys:
            return max(keys)
        return -1

    def get_keys(self, rev=False):
        keys = list(self.container.keys())
        if rev:
            keys = keys[::-1]
        keys.sort()
        return keys

    def __getitem__(self, i):
        if i == -1:
            i = self.get_max_degree()
        if i not in self.container:
            return 0
        return self.container[i]

    def __setitem__(self, i, a):
        if i < 0:
            raise FPNegativeDegreeNotAllowed()
        if a == 0:
            if i in self.container:
                del self.container[i]
        else:
            self.container[i] = a

    def _element_str(self, idx, empty):
        element = self.container[idx]
        if element:
            return "{}{}x^{} ".format(
                get_sign(element, empty), abs(element) if element != 1 else "", idx
            )
        return ""

    def __str__(self):
        s = ""
        keys = self.get_keys()
        if keys:  # Empty polynom should return 0
            if len(keys) != 1:
                # Largest degree
                e = self.container[keys[-1]]
                if e:
                    s = self._element_str(keys[-1], True)

                # Intermediate
                for i in keys[1:-1][::-1]:
                    s += self._element_str(i, not s)

            # Lowest degree
            if keys[0] == 0:
                e = self.container[0]
                if e or not s:
                    s += "{}{}".format(get_sign(e, not s), abs(e))
            else:
                s += self._element_str(keys[0], not s)
            return s
        else:
            return "0"


if __name__ == "__main__":
    p = FastPolynom()
    p[1] = 1
    p[300] = 20
    p[34] = 4
    p[12] = 8
    p[55] = 5
    print(p[-1])
    print(p)
    print(type(p) == FastPolynom)

# -> max_degree(n), O(1), O(1)
# -> n, O(1), O(n)
# -> n, O(log(n)), O(1)
