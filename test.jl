data = Array{Float32}(undef,n)
read!("/Users/joshuaighalo/Downloads/raw/0_1_12072018_1206/0_1_12072018_1206.bin", data)



test_2 = (read("/Users/joshuaighalo/Downloads/raw/0_1_12072018_1206/0_1_12072018_1206.bin"))
len_test_2 = length(test_2)
dataCols = 11
dataRows = Int(length(test_2)/dataCols)
coll_data = reshape(test_2, (dataRows, dataCols))
