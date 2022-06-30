// Magma code

// Returns an integer sequence as a string
function IntegerSequenceToString(s)
	string:="[";
	for i in [1..#s] do
		string:=string cat IntegerToString(s[i]);
		if i lt #s then
			string:=string cat ",";
		else
			string:=string cat "]";
		end if;
	end for;
	return string;
end function;

// Returns a sequence of vertices as a string
function VerticesToString(s)
	string:="[";
	for i in [1..#s] do
		string:=string cat IntegerSequenceToString(s[i]);
		if i lt #s then
			string:=string cat ",";
		else
			string:=string cat "]";
		end if;
	end for;
	return string;
end function;

D:=3;

A:=[
	[[0,0,0],[1,0,0],[1,2,0],[2,2,0],[1,0,2],[2,0,2],[2,2,2],[3,2,2]],
	[[0,0,0],[2,0,0],[2,4,0],[2,0,4]],
	[[0,0,0],[1,0,0],[2,4,0],[3,0,4]],
	[[0,0,0],[1,0,0],[2,3,0],[1,0,3],[-1,-6,3]],
	[[0,0,0],[3,0,0],[0,3,0],[0,0,3]],
	[[0,0,0],[1,0,0],[1,4,0],[3,0,4],[-1,4,-4]],
	[[0,0,0],[1,0,0],[4,6,0],[4,0,6]],
	[[0,0,0],[1,0,0],[2,3,0],[5,3,9]],
	[[0,0,0],[1,0,0],[3,4,0],[7,4,8]],
	[[0,0,0],[1,0,0],[2,3,0],[1,0,3],[2,0,3],[3,3,3]],
	[[0,0,0],[1,0,0],[2,5,0],[3,0,5]],
	[[0,0,0],[1,0,0],[1,2,0],[3,2,4],[2,2,0],[4,2,4]]
];

N:=Max([NumberOfPoints(Polytope(s)) : s in A]);

list:=[[] : i in [1..N]];

for s in A do
	P:=Polytope(s);
	n:=NumberOfPoints(P);
	Append(~list[n],Polytope([Eltseq(p) : p in AffineNormalForm(P)]));
end for;

for i in [N..1 by -1] do
	if #list[i] eq 0 then break; end if;
	for P in list[i] do
		for v in Vertices(P) do
			Q:=Polytope(Exclude(Points(P),v));
			if Dimension(Q) eq D and Width(Q) gt 1 then
				QQ:=Polytope([Eltseq(p) : p in AffineNormalForm(Q)]);
				Append(~list[NumberOfPoints(QQ)],QQ);
			end if;
		end for;
	end for;
	set_p:=SequenceToSet(list[i-1]);
	list[i-1]:=SetToSequence(set_p);
	printf "%o - %o\n",i,#list[i];
end for;

set:=SequenceToSet(&cat(list));
#set;

function l2(P)
    n1:=#InteriorPoints(2*P);
    n2:=0;
    for F in Facets(P) do
        n2:=n2+#InteriorPoints(F);
    end for;
    return n1-n2;
end function;

for P in set do
    if l2(P) eq 0 then
        print "Exception!";
    end if;
end for;
