# Seasonality data processing in FinRegistry
# Modified from endpoints_monthly_counts_FINNGEN.R

library("data.table")
library("lubridate")

# Read birth dates

print("Reading birth date data")
path_birth_dates = "/data/processed_data/minimal_phenotype/minimal_phenotype_main_dec2022.csv"
birth_dat = fread(path_birth_dates, select = c("FINREGISTRYID", "date_of_birth"))

# Reformat dates

print("Formatting birth dates")
birth_dates = dmy(birth_dat$date_of_birth)

# Fix invalid dates
# All are 29.2. for years that are not leap years
# --> replace with 28.2.

invalid_dates = birth_dat$date_of_birth[is.na(birth_dates)]
substring(invalid_dates, 1, 2) = "28"
fixed_dates = dmy(invalid_dates)
birth_dates[is.na(birth_dates)] = fixed_dates

# Add FinRegistry ID to birth dates

names(birth_dates) = birth_dat$FINREGISTRYID

# Process longitudinal endpoints data
# OMITed endpoints not included

print("Setting up the connection")
path_endpoint = "/data/processed_data/endpointer/R10/longitudinal_endpoints_no_omits_DF10_2022_09_29.txt.ALL.gz"
con_endpoint = gzcon(file(path_endpoint, "rb"))
col_names = unlist(strsplit(readLines(con_endpoint, n = 1), split = ","))

first_event_summary = setDT(NULL)
last_iter_data = data.table(FINREGISTRYID=character(), ENDPOINT=character())
batch_size = 10e6 # approx 486 iterations in FinRegistry data

time = Sys.time()
iteration = 1

while(TRUE){
  print(paste("Iteration:", iteration))
  
  # Read chunk
  dat_endpoint_batch = fread(text=readLines(con_endpoint, n = batch_size), sep = ",")
  
  # Break loop if all data has been read
  if(nrow(dat_endpoint_batch) == 0) break
  
  # Unique extracts first pair (FINREGISTYID, ENDPOINT) as table is sorted w.r.t. time
  names(dat_endpoint_batch) = col_names
  first_event_batch = unique(dat_endpoint_batch, by = c("FINREGISTRYID", "ENDPOINT"))
  
  # Add event date
  event_date = ymd(birth_dates[first_event_batch$FINREGISTRYID]) + dyears(first_event_batch$EVENT_AGE)
  first_event_batch$EVENT_MONTH = month(event_date)
  first_event_batch$EVENT_YEAR = year(event_date)
  
  # Anti-join with last individual in the previous iteration to avoid counting (FINREGISTRYID, event) pairs
  first_event_batch = first_event_batch[!last_iter_data, on = c("FINREGISTRYID", "ENDPOINT")]
  
  # Store data on last individual to avoid double counting in next iteration
  last_iter_data = first_event_batch[FINREGISTRYID == tail(FINREGISTRYID, 1), ] # error here?
  
  # Aggregate data
  first_event_summary_batch = first_event_batch[, .(COUNT=.N), by = c("ENDPOINT", "EVENT_YEAR", "EVENT_MONTH")]
  
  # Combine data
  first_event_summary = rbind(first_event_summary, first_event_summary_batch)
  first_event_summary = first_event_summary[, .(COUNT=sum(COUNT)), by = c("ENDPOINT", "EVENT_YEAR", "EVENT_MONTH")]
  
  iteration = iteration + 1
}

# Close connection
close(con_endpoint)

print("Time passed:")
print(Sys.time() - time)

# Remove personal data 
df_red = first_event_summary
df_green = first_event_summary[first_event_summary$COUNT >= 5, ]

print("Writing the output")
output_path_red = "/data/projects/seasonality/results/FINREGISTRY_endpoints_monthly_count_red.txt"
output_path_green = "/data/projects/seasonality/results/FINREGISTRY_endpoints_monthly_count_green.txt"

fwrite(df_red, file = output_path_red, sep = "\t", quote = F)
fwrite(df_green, file = output_path_green, sep = "\t", quote = F)
