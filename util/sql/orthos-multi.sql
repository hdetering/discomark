select count(*) 
from (
  select id_ortholog, count(distinct id_species) as c 
  from sequences 
  group by id_ortholog having c>1
)