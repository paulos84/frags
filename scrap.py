with open('C:\\Users\\Paul\\M1_book.txt', encoding="utf8") as w:
    r = w.read()

import re

m = re.findall(r'\[[^\[\]]*\d+(?:-\d+)+\]', r)
lst = []
for a in m:
    lst.append(a)

print(lst)
print(len(lst))
#     lst.append(re.findall(r'\d+(?:-\d+)+',a))
# [a[0] for a in lst if len(a) >0]
# cases = [a[0] for a in lst if len(a) >0]


# m = re.findall('\[[^\[\]]*\]', r)

"""
m = re.findall(r'\[[^\[\]]*\d+(?:-\d+)+\]', r)
lst = []
for a in m:
    lst.append(a)

print(lst)
print(len(lst))
"""


