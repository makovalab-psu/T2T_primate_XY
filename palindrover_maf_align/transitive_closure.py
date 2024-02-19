#!/usr/bin/env python3
"""
Cluster items by transitive closure of pairwise associations.

Conceptually, we are given a series of pairwise "connections" -- a pair of
items which should be in the same cluster. We group any items for which there
is a pairwise connection in a cluster.

Note: this is a quick-and-dirty implementation, and could be improved. We
maintain a directed graph (represented by a dict) which is mainly cycle-free.
The only cycles have length 1; these are called roots. Each connected component
in the true cluster graph is represented by a tree in which all paths lead to
the root. The "connect" operation is implemented by finding the root of each of
the two items and changing one root to point to the other. While this works
well for small graphs, larger graphs will undoubtedly expose issues w.r.t.
runtime.
"""

from collections import defaultdict


class TransitiveClosure(object):

	def __init__(self):
		self.itemTowardRoot = {}
		self.mark_as_changed()

	def mark_as_changed(self):
		self.itemToGroup = None
		self.clusterIds = None

	def add(self,item):
		if (item not in self.itemTowardRoot[item]):
			self.itemTowardRoot[item] = item
		self.mark_as_changed()

	def connect(self,item1,item2):
		root1 = self.root_of(item1)
		root2 = self.root_of(item2)
		if (root1 == None): root1 = self.itemTowardRoot[item1] = item1
		if (root2 == None): root2 = self.itemTowardRoot[item2] = item2
		self.itemTowardRoot[root2] = root1
		self.mark_as_changed()

	def derive_clusters(self):
		if (self.itemToGroup != None): return self.itemToGroup
		self.clusterIds = None

		rootToItems = defaultdict(list)
		for item in self.itemTowardRoot:
			root = self.root_of(item)
			rootToItems[root] += [item]

		self.itemToGroup = {}
		for root in rootToItems:
			minItem = min(rootToItems[root])
			self.itemToGroup[minItem] = rootToItems[root]

		return self.itemToGroup

	def ordered_cluster_ids(self): # ordered by decreasing size
		if (self.itemToGroup == None): self.derive_clusters()

		items = [(-len(self.itemToGroup[item]),item) for item in self.itemToGroup]
		items.sort()
		self.clusterIds = [item for (_,item) in items]

		return self.clusterIds

	def root_of(self,item):
		# (this is intended to be a private method)
		if (item not in self.itemTowardRoot):
			return None

		(prevItem,nextItem) = (item,self.itemTowardRoot[item])
		while (nextItem != prevItem):
			(prevItem,nextItem) = (nextItem,self.itemTowardRoot[nextItem])

		return nextItem

