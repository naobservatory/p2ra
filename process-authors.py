#!/usr/bin/env python3

import re
from collections import defaultdict

by_last_name = defaultdict(set)
all_authors = set()
# Manually copied from papers
with open("authors.txt") as inf:
    for line in inf:
        line = line.removesuffix("\n")
        paper, authors = line.split(": ")

        for author in authors.split(","):
            author = author.strip()
            author = author.removeprefix("and ")

            for f, r in [
                    (".", ""),
                    ("Thomas Nordahl Petersen", "Thomas N Petersen"),
                    ("Christina Aaby Svendsen", "Christina A Svendsen"),
                    ]:
                author = author.replace(f,r)

            # remove middle names
            while re.match(".* [A-Z] .*", author):
                author = re.sub("(.+) [A-Z] (.+)", r"\1 \2", author)

            by_last_name[author.split(" ")[-1]].add(author)
            all_authors.add(author)

for last_name in sorted(by_last_name):
    if len(by_last_name[last_name]) == 1: continue
    print(last_name)
    for author in sorted(by_last_name[last_name]):
        print(" ", author)

print(", ".join(sorted(all_authors)))
