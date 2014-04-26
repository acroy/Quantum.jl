

s = d.State("sr")
a = d.Basis("a", [1:10])
b = d.Basis("b", ["$i" for i=1:4]);
c = d.Basis("c", [1:20]);
r = d.tensor(a,b,c)
sr = d.normalize!(d.StateRep(s, [1 for i=1:length(r)], r))
op = (sr*sr')


q = d.Basis("q", b[1:2])
sq = d.StateRep(s, [1, 1], q)
sq = sq*sq
sq[1] = 0
sq[4] = 0
d.normalize!(sq)

print("")