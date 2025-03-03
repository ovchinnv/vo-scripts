# Remove null entries from alt_id because matlab cannot tolerate it
#
.open /home/taly/flurepo/flu_strains.db
.tables

update flu_strains90 set msa="" where msa is NULL ;
update flu_strains95 set msa="" where msa is NULL ;
update flu_strains97 set msa="" where msa is NULL ;
