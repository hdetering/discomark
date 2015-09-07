select count(*) 
from (
  select s.id_ortholog, count(distinct id_species) as c 
  from 
    sequences s inner join 
    primer_sets ps on ps.id_ortholog = s.id_ortholog 
  group by s.id_ortholog 
  having c>1
)