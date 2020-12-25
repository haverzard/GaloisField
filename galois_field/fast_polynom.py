#!/usr/bin/python3

from .exceptions import FPNegativeDegreeNotAllowed
from .util import get_sign


class FastPolynom:
    """
    Specialized polynom for galois field
    """

    def __init__(self, container=None):
        self.cache = False
        self.cached_keys = None
        self.container = container or {}

        # Cleanup zeros
        if container:
            keys = list(container.keys())
            for x in keys:
                if not container[x]:
                    del container[x]

    def deepcopy(self):
        return FastPolynom(self.container.copy())

    def broadcast_modulo(self, p):
        for i in self.container:
            self.container[i] %= p

    def get_max_degree(self):
        # Check for cache
        if self.cache:
            return self.cached_keys[-1]

        # Get max degree
        keys = list(self.container.keys())
        if keys:
            return max(keys)
        return -1

    def get_keys(self, rev=False):
        # Check for cache
        if self.cache:
            return self.cached_keys

        # Get keys
        keys = list(self.container.keys())
        if rev:
            keys = keys[::-1]
        keys.sort()

        # Set cache
        self.cached_keys = keys
        self.cache = True

        return keys

    def compute(self, x):
        total = 0
        for i in self.get_keys(rev=True):
            total = total * x + self[i]
        return total

    def __getitem__(self, i):
        # Special case: get max degree's element
        if i == -1:
            i = self.get_max_degree()

        # Get element from degree i
        if i not in self.container:
            return 0
        return self.container[i]

    def __setitem__(self, i, a):
        if i < 0:
            raise FPNegativeDegreeNotAllowed()
        if a == 0:
            if i in self.container:
                self.cache = False
                del self.container[i]
        else:
            self.cache = self.cache and i in self.container
            self.container[i] = a

    def _element_str(self, idx, empty):
        element = self.container[idx]
        if element:
            return "{}{}x^{} ".format(
                get_sign(element, empty), abs(element) if element != 1 else "", idx
            )
        return ""

    def __eq__(self, x):
        if self.get_keys() != x.get_keys():
            return False
        for i in self.get_keys():
            if self[i] != x[i]:
                return False
        return True

    def __ne__(self, x):
        return not (self == x)

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
