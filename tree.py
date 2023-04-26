from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Generic, TypeVar

T = TypeVar("T")
S = TypeVar("S")


@dataclass
class Tree(Generic[T]):
    data: T
    children: list[Tree[T]] = field(default_factory=list)

    def _str_helper(self, depth: int) -> str:
        _spacer = "."
        return f"{_spacer * depth}{self.data}\n" + "".join(
            c._str_helper(depth + 1) for c in self.children
        )

    def __str__(self) -> str:
        return self._str_helper(0)

    def __iter__(self):
        yield self
        for child in self.children:
            yield from child

    def __getitem__(self, val: T) -> Tree[T] | None:
        return self._get_subtree(val)

    def __contains__(self, item: T) -> bool:
        return not (self[item] is None)

    def _get_subtree(self, val: T) -> Tree[T] | None:
        """Depth-first search for taxid"""
        for subtree in self:
            if subtree.data == val:
                return subtree
        else:
            return None

    def to_list(self) -> list:
        return [self.data] + [c.to_list() for c in self.children]

    def map(self, f: Callable[[T], S]) -> Tree[S]:
        return Tree(f(self.data), [c.map(f) for c in self.children])

    @staticmethod
    def tree_from_list(input: list) -> Tree:
        return Tree(
            data=input[0], children=[Tree.tree_from_list(c) for c in input[1:]]
        )
