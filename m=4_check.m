groups := ExtensionsOfElementaryAbelianGroup(2,4,PGL(4,2));
for G in groups do
    table := CharacterTable(PermutationGroup(G));
    codegrees := [1];
    i := 2;
    orderG := Order(G);
    while i le #table do
        cod := orderG  / Order(Kernel(table[i])) / table[i][1];
        codegrees := Append(codegrees, cod);
        i := i + 1;
    end while;
    codset := [];
    for j in Sort(codegrees) do
        if not j in codset then
            codset := Append(codset, j);
        end if;
    end for;
    print codset;
end for;