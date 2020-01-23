# Queries for the database

source("r1_db_funcs.R")

library(plyr)
library(DBI)
library(RSQLite)

drv <- dbDriver("SQLite")
data <- dbConnect(drv,"tsetse2018.db") # connect to the database

treatment <- query.func(query="SELECT * from treatment")

#******data on adult deaths**************
adult.deaths <- query.func(query="SELECT T.name, T.starting_size, T.date_of_emergence_min, A.adults_id, A.date_of_death, A.weekend_death, A.tray, A.colony_tray
                           FROM treatment T, adults A
                           WHERE T.treatment_id=A.treatment_id")

#*****linking mother and larvipositions****
mother.larvipositions <- query.func(query="SELECT T.name,T.date_of_emergence_min,  A.adults_id, A.well_number, A.date_of_death, A.colony_tray,A.hatchet_length, L.larviposition_id, L.larviposition_date, L.larviposition_number, L.weekend, L.wet_weight, L.dry_weight, L.residual_dry_weight, L.abortion_stage, L.abortion
                             FROM treatment T, adults A, larviposition L
                             WHERE T.treatment_id=A.treatment_id AND A.adults_id=L.adults_id")

#*****pupae******
pupae <- query.func(query="SELECT T.name, T.date_of_emergence_min, A.adults_id, A.well_number, L.larviposition_id, L.larviposition_date, L.larviposition_number, L.weekend, L.abortion_stage, L.abortion
                             ,L.wet_weight,P.pupa_length,P.pupa_width,P.emerged,P.date_emerged,P.date_of_death,P.weekend_death,P.sex,P.treatment, P.weight_emerged
                             FROM treatment T, adults A, larviposition L, offspring P
                             WHERE T.treatment_id=A.treatment_id AND A.adults_id=L.adults_id AND L.larviposition_id=P.larviposition_id")


dbDisconnect(data)
