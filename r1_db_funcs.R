library("RSQLite") 

#**********************FUNCTION FOR RETRIEVING DATA FROM DATABASE********************************
# function used when importing data to first retrieve the last primary
# keys for each table so that the next set of primary key numbers are sequential 

query.func <- function(db=data,query){ # function for retrieving data from the database
  sendquery <- dbSendQuery(db,query)
  sendquery
  table <- fetch(sendquery,n=-1)
  dbClearResult(sendquery)
  return(table)
} 

#***********************FUNCTION FOR GENERATING PRIMARY KEYS*********************************
# arguments- query =  retrieve the orimary key used in that table, 
# col = the number of entries to be added
primary_key.func<- function(query="SELECT * FROM X",
                            col=gps_table$date){
  current_id <- as.numeric(query.func(db=data,query)) 
  current_id[current_id %in% NA] <- 0
  next_id <- current_id + 1 
  new_id <- seq(next_id,next_id + length(col)-1,1) 
  return(new_id)
} 
#################################################################################################