# ----------------------------------------------------------------------------------- #
# Copyright (c) 2016 Varnerlab
# Robert Frederick School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #

# include -
#time_start = 0.0
#time_stop = 120.0
#time_step_size = 0.01
#path_to_par = "randomized_parameter_ensemble.dat"
#path_to_error = "error_archive_random_ensemble.dat"
#path_to_simulation = "./simulation/"
#number_of_samples = 10

function strictly_sample_ensemble(time_start,time_stop,time_step_size,path_to_par,path_to_error,path_to_simulation,number_of_samples)

    # check file ok?
    is_file_path_ok(path_to_par)
    is_file_path_ok(path_to_error)
    is_dir_path_ok(path_to_simulation)
    # check time ok?
    is_time_ok(time_start)
    is_time_ok(time_stop)
    is_time_ok(time_step_size)
    is_number_ok(number_of_samples)
    is_time2_greater_than_time1(time_start,time_stop)

    # Script to solve the balance equations -
    # Load the data dictionary -
    data_dictionary = DataDictionary(time_start,time_stop,time_step_size)

    # get my initial parameter guess from the previous run -
    par_array = readdlm(path_to_par)
    error_archive = vec(readdlm(path_to_error))
    idx_error_sort = sortperm(error_archive,rev=false)
    is_number2_bigger_than_number1(number_of_samples,size(par_array)[2])
    par_array = par_array[:,idx_error_sort[1:number_of_samples]]
    (number_of_rows,number_of_samples) = size(par_array)

    # Declare a progress meter for user feedback -
    p = Progress(number_of_samples,color=:yellow)

    # main loop -
    for outer_sample_index = 1:number_of_samples
        # what parameters are we looking at?
        parameter_set = par_array[:,outer_sample_index]
        # make a copy of the data dictionary -
        copy_data_dictionary = deepcopy(data_dictionary)

        # how many different types of parameters do we have?
        parameter_mapping_array = copy_data_dictionary["parameter_name_mapping_array"]
        number_of_binding_parameters = length(copy_data_dictionary["binding_parameter_dictionary"])
        number_of_control_parameters = length(copy_data_dictionary["control_parameter_dictionary"])

        # Update data dictionary to match new parameters before running the simulation -
        for index = 1:length(parameter_mapping_array)

            parameter_name = parameter_mapping_array[index]

            if index <= number_of_binding_parameters
                copy_data_dictionary["binding_parameter_dictionary"][parameter_name] = parameter_set[index]
            elseif (index>number_of_binding_parameters && index<=(number_of_binding_parameters+number_of_control_parameters))
                copy_data_dictionary["control_parameter_dictionary"][parameter_name] = parameter_set[index]
            else
                copy_data_dictionary[parameter_name] = parameter_set[index]
            end
        end

        # run the simulation -
        (time_array,simulation_state_array) = add_atra_simulation(time_start,time_stop,time_step_size,copy_data_dictionary);
        local_data = [time_array simulation_state_array];
        data_filename = path_to_simulation*"sim_data_"*string(outer_sample_index)*".dat"

        #data_filename = "./simulation/sim_data_"*string(outer_sample_index)*".dat"
        writedlm(data_filename,local_data);

        message_string = "Completed $(outer_sample_index) of $(number_of_samples)"
        # update the progress bar -
        ProgressMeter.next!(p; showvalues = [(:status,message_string)])
        #println(message_string)
    end
end
