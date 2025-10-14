# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

__name__ = "biotite.sequence.phylo"
__author__ = "Patrick Kunzmann, Tom David Müller"
__all__ = ["Tree", "TreeNode", "as_binary", "TreeError"]

import copy as copy_module
import numpy as np

from typing import List, Optional


class Tree:
    """A rooted tree representation"""
    
    _root: object
    _leaves: List[object]
    
    def __init__(self, root):
        if root is None:
            raise ValueError("Root cannot be None")
        root.as_root()
        self._root = root
        
        leaves_unsorted = self._root.get_leaves()
        leaf_count = len(leaves_unsorted)
        indices = np.array([leaf.index for leaf in leaves_unsorted])
        self._leaves = [None] * leaf_count
        
        for i in range(len(indices)):
            index = int(indices[i])
            if index >= leaf_count or index < 0:
                raise TreeError("The tree's indices are out of range")
            self._leaves[index] = leaves_unsorted[i]
    
    @property
    def root(self):
        return self._root
    
    @property
    def leaves(self):
        return copy_module.copy(self._leaves)
    
    def as_graph(self):
        """Obtain a graph representation of the Tree using NetworkX"""
        try:
            import networkx as nx
        except ImportError:
            raise ImportError("networkx is required for as_graph()")
        
        graph = nx.DiGraph()
        node_repr: dict = {}

        queue: List[object] = copy_module.copy(self._leaves)
        queue_set: set = set(self._leaves)
        
        while len(queue) > 0:
            node = queue.pop(0)
            
            if node.is_leaf():
                node_repr[node] = node.index
            else:
                children = node.children
                children_handled = True
                for child in children:
                    if child not in node_repr:
                        children_handled = False
                
                if not children_handled:
                    queue.append(node)
                    continue
                else:
                    repr_val = tuple(node_repr[child] for child in children)
                    node_repr[node] = repr_val
                    for child in children:
                        graph.add_edge(
                            repr_val, node_repr[child], distance=child.distance
                        )
            
            if not node.is_root():
                parent = node.parent
                if parent not in queue_set:
                    queue.append(parent)
                    queue_set.add(parent)
            
            queue_set.remove(node)
        
        return graph

    def get_distance(self, index1: int, index2: int, topological: bool = False):
        """Get the distance between two leaf nodes"""
        return self._leaves[index1].distance_to(
            self._leaves[index2], topological
        )
        if topological:
            return int(distance_in_edges)  # convert to int
        return float(distance_in_branch_lengths)
    
    def to_newick(self, labels = None, include_distance: bool = True, 
                  round_distance = None):
        """Obtain the Newick notation of the tree"""
        root_obj = self._root
        return root_obj.to_newick(
            labels, include_distance, round_distance
        ) + ";"
    
    @staticmethod
    def from_newick(newick: str, labels = None):
        """Create a tree from a Newick notation"""
        newick = newick.strip()
        if len(newick) == 0:
            raise ValueError("Newick string is empty")
        if newick[-1] == ";":
            newick = newick[:-1]
        root, distance = TreeNode.from_newick(newick, labels)
        return Tree(root)

    def __str__(self):
        return self.to_newick()
    
    def __len__(self):
        return len(self._leaves)
    
    def __eq__(self, item):
        if not isinstance(item, Tree):
            return False
        return self._root == item._root
    
    def __hash__(self):
        return hash(self._root)


class TreeNode:
    """A node in a rooted tree"""

    _index: int
    _distance: float
    _is_root: bool
    _parent: object
    _children: object

    def __init__(self, children = None, distances = None, index = None):
        self._is_root = False
        self._distance = 0.0
        self._parent = None
        self._children = None
        self._index = -1
        
        if index is None:
            # Node is intermediate -> has children
            if children is None or distances is None:
                raise TypeError(
                    "Either reference index (for terminal node) or "
                    "child nodes including the distance "
                    "(for intermediate node) must be set"
                )
            for item in children:
                if not isinstance(item, TreeNode):
                    raise TypeError(
                        f"Expected 'TreeNode', but got '{type(item).__name__}'"
                    )
            for item in distances:
                #if not isinstance(item, float) and not isinstance(item, int):
                #if not isinstance(item, float):
                if not isinstance(item, (float, int)):
                    raise TypeError(
                        f"Expected 'float' or 'int', "
                        f"but got '{type(item).__name__}'"
                    )
            if len(children) == 0:
                raise TreeError(
                    "Intermediate nodes must at least contain one child node"
                )
            if len(children) != len(distances):
                raise ValueError(
                    "The number of children must equal the number of distances"
                )
            for i in range(len(children)):
                for j in range(len(children)):
                    if i != j and children[i] is children[j]:
                        raise TreeError(
                            "Two child nodes cannot be the same object"
                        )
            self._index = -1
            self._children = tuple(children)
            for child, distance in zip(children, distances):
                child._set_parent(self, float(distance))
        elif index is not None and index < 0:
            raise ValueError("Index cannot be negative")
        else:
            # Node is terminal -> has no children
            if children is not None or distances is not None:
                raise TypeError(
                    "Reference index and child nodes are mutually exclusive"
                )
            self._index = index if index is not None else -1
            self._children = None
    
    def _set_parent(self, parent, distance: float):
        if parent is None:
            raise ValueError("Parent cannot be None")
        if self._parent is not None or self._is_root:
            raise TreeError("Node already has a parent")
        self._parent = parent
        self._distance = distance
    
    def copy(self):
        """Create a deep copy of this TreeNode"""
        if self.is_leaf():
            return TreeNode(index=self._index)
        else:
            distances = [child.distance for child in self._children]
            children_clones = [child.copy() for child in self._children]
            return TreeNode(children_clones, distances)

    @property
    def index(self):
        return None if self._index == -1 else self._index
    
    @property
    def children(self):
        return self._children
    
    @property
    def parent(self):
        return self._parent
    
    @property
    def distance(self):
        return None if self._parent is None else self._distance

    def is_leaf(self):
        """Check if the node is a leaf node"""
        return False if self._index == -1 else True
    
    def is_root(self):
        """Check if the node is a root node"""
        return bool(self._is_root)
    
    def as_root(self):
        """Convert the node into a root node"""
        if self._parent is not None:
            raise TreeError("Node has parent, cannot be a root node")
        self._is_root = True
    
    def distance_to(self, node, topological: bool = False):
        """Get the distance of this node to another node"""
        distance = 0.0
        lca = self.lowest_common_ancestor(node)
        if lca is None:
            raise TreeError("The nodes do not have a common ancestor")
        current_node = self
        while current_node is not lca:
            if topological:
                distance += 1.0
            else:
                distance += current_node._distance
            current_node = current_node._parent
        current_node = node
        while current_node is not lca:
            if topological:
                distance += 1.0
            else:
                distance += current_node._distance
            current_node = current_node._parent
        return distance
    
    def lowest_common_ancestor(self, node):
        """Get the lowest common ancestor of this node and another node"""
        lca = None
        self_path = _create_path_to_root(self)
        other_path = _create_path_to_root(node)
        for i in range(-1, -min(len(self_path), len(other_path))-1, -1):
            if self_path[i] is other_path[i]:
                lca = self_path[i]
            else:
                break
        return lca
    
    def get_indices(self):
        """Get an array of reference indices for leaf nodes"""
        return np.array(
            [leaf._index for leaf in self.get_leaves()], dtype=np.float64
        )

    def get_leaves(self):
        """Get a list of leaf nodes that are child nodes of this node"""
        leaf_list: List[object] = []
        _get_leaves(self, leaf_list)
        return leaf_list
    
    def get_leaf_count(self):
        """Get the number of direct or indirect leaves of this node"""
        return _get_leaf_count(self)
    
    def to_newick(self, labels = None, include_distance: bool = True, 
                  round_distance = None):
        """Obtain the node represented in Newick notation"""
        if self.is_leaf():
            if labels is not None:
                label = labels[self._index]
                illegal_chars = [",", ":", ";", "(", ")"]
                for char in illegal_chars:
                    if char in label:
                        raise ValueError(
                            f"Label '{label}' contains "
                            f"illegal character '{char}'"
                        )
            else:
                label = str(self._index)
            if include_distance:
                if round_distance is None:
                    return f"{label}:{self._distance}"
                else:
                    return f"{label}:{self._distance:.{round_distance}f}"
            else:
                return f"{label}"
        else:
            child_strings = [child.to_newick(
                labels, include_distance, round_distance
            ) for child in self._children]
            if include_distance:
                if round_distance is None:
                    return f"({','.join(child_strings)}):{self._distance}"
                else:
                    return (
                        f"({','.join(child_strings)}):"
                        f"{self._distance:.{round_distance}f}"
                    )
            else:
                return f"({','.join(child_strings)})"
    
    @staticmethod
    def from_newick(newick: str, labels = None):
        """Create a node and all its child nodes from a Newick notation"""
        newick = "".join(newick.split())

        subnewick_start_i = -1
        subnewick_stop_i = -1
        for i in range(len(newick)):
            char = newick[i]
            if char == "(":
                subnewick_start_i = i
                break
            if char == ")":
                raise ValueError("Bracket closed before it was opened")
        
        for i in reversed(range(len(newick))):
            char = newick[i]
            if char == ")":
                subnewick_stop_i = i + 1
                break
            if char == "(":
                raise ValueError("Bracket was opened but not closed")
        
        if subnewick_start_i == -1 and subnewick_stop_i == -1:
            # No brackets -> Leaf node
            label_and_distance = newick
            try:
                label, distance_str = label_and_distance.split(":")
                distance = float(distance_str)
            except ValueError:
                distance = 0.0
                label = label_and_distance
            index = int(label) if labels is None else labels.index(label)
            return TreeNode(index=index), distance
        
        else:
            # Intermediate node
            if subnewick_stop_i == len(newick):
                label = None
                distance = 0.0
            else:
                label_and_distance = newick[subnewick_stop_i:]
                try:
                    label, distance_str = label_and_distance.split(":")
                    distance = float(distance_str)
                except ValueError:
                    distance = 0.0
                    label = label_and_distance
                distance = float(distance)
            
            subnewick = newick[subnewick_start_i + 1 : subnewick_stop_i - 1]
            if len(subnewick) == 0:
                raise ValueError(
                    "Intermediate node must at least have one child"
                )
            
            comma_pos: List[int] = []
            level = 0
            for i, char in enumerate(subnewick):
                if char == "(":
                    level += 1
                elif char == ")":
                    level -= 1
                elif char == ",":
                    if level == 0:
                        comma_pos.append(i)
                if level < 0:
                    raise ValueError(
                        "Bracket closed before it was opened"
                    )
        
            children: List[object] = []
            distances: List[float] = []
            for i, pos in enumerate(comma_pos):
                if i == 0:
                    child, dist = TreeNode.from_newick(
                        subnewick[:pos], labels=labels
                    )
                else:
                    prev_pos = comma_pos[i - 1]
                    child, dist = TreeNode.from_newick(
                        subnewick[prev_pos + 1 : pos], labels=labels
                    )
                children.append(child)
                distances.append(dist)
            
            if len(comma_pos) != 0:
                child, dist = TreeNode.from_newick(
                    subnewick[comma_pos[-1] + 1:], labels=labels
                )
            else:
                child, dist = TreeNode.from_newick(
                    subnewick, labels=labels
                )
            children.append(child)
            distances.append(dist)
            return TreeNode(children, distances), distance

    def __str__(self):
        return self.to_newick()
    
    def __eq__(self, item):
        if not isinstance(item, TreeNode):
            return False
        node = item
        if self._distance != node._distance:
            return False
        if self._index != -1:
            if self._index != node._index:
                return False
        else:
            if frozenset(self._children) != frozenset(node._children):
                return False
        return True
    
    def __hash__(self):
        children_set = frozenset(self._children) \
                       if self._children is not None else None
        return hash((self._index, children_set, self._distance))
    

    def __eq__(self, other):
        if not isinstance(other, TreeNode):
            return False
        if getattr(self, "index", None) != getattr(other, "index", None):
            return False
        if len(self._children) != len(other._children):
            return False
        for c1, c2 in zip(self._children, other._children):
            if c1 != c2:
                return False
        if hasattr(self, "_distances") and hasattr(other, "_distances"):
            if any(abs(d1 - d2) > 1e-6 for d1, d2 in zip(self._distances, other._distances)):
                return False
        return True


def _get_leaves(node, leaf_list: List[object]):
    if node._index == -1:
        for child in node._children:
            _get_leaves(child, leaf_list)
    else:
        leaf_list.append(node)


def _get_leaf_count(node):
    count = 0
    if node._index == -1:
        for child in node._children:
            count += _get_leaf_count(child)
        return count
    else:
        return 1


def _create_path_to_root(node):
    """Create a list of nodes representing the path from this node to the root"""
    path: List[object] = []
    current_node = node
    while current_node is not None:
        path.append(current_node)
        current_node = current_node._parent
    return path


def as_binary(tree_or_node):
    """Convert a tree into a binary tree"""
    if isinstance(tree_or_node, Tree):
        node, _ = _as_binary(tree_or_node.root)
        return Tree(node)
    elif isinstance(tree_or_node, TreeNode):
        node, _ = _as_binary(tree_or_node)
        return node
    else:
        raise TypeError(
            f"Expected 'Tree' or 'TreeNode', not {type(tree_or_node).__name__}"
        )


def _as_binary(node):
    """The actual logic for as_binary()"""
    children = node.children
    if children is None:
        return TreeNode(index=node.index), node.distance
    elif len(children) == 1:
        child, distance = _as_binary(node.children[0])
        if node.is_root():
            return child, None
        else:
            child_dist = distance if distance is not None else 0.0
            return child, node.distance + child_dist
    elif len(children) > 2:
        results = [_as_binary(child) for child in children]
        rem_children: List[object] = [result[0] for result in results]
        distances_list: List[float] = [result[1] if result[1] is not None else 0.0 for result in results]
        current_div_node = None
        
        while len(rem_children) > 0:
            if current_div_node is None:
                current_div_node = TreeNode(
                    rem_children[:2],
                    distances_list[:2]
                )
                rem_children.pop(0)
                rem_children.pop(0)
                distances_list.pop(0)
                distances_list.pop(0)
            else:
                current_div_node = TreeNode(
                    [current_div_node, rem_children[0]],
                    [0.0, distances_list[0]] 
                )
                rem_children.pop(0)
                distances_list.pop(0)
        return current_div_node, node.distance
    else:
        results = [_as_binary(child) for child in children]
        binary_children: List[object] = [result[0] for result in results]
        distances_list: List[float] = [result[1] if result[1] is not None else 0.0 for result in results]
        return TreeNode(binary_children, distances_list), node.distance


class TreeError:
    """An exception that occurs in context of tree topology"""
    message: str
    
    def __init__(self, message: str):
        self.message = message
    
    def __str__(self):
        return self.message