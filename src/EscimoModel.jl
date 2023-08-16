  using ModelingToolkit
  using DifferentialEquations
  import ModelingToolkit.IfElse
  import OMRuntimeExternalC
  include("utilityFunctions.jl")
  Symbolics.@register_symbolic ESCIMO_IF_THEN_ELSE(condition, result_true, result_false)
  Symbolics.@register_symbolic ESCIMO_Population_Lookup_bn(y)
  Symbolics.@register_symbolic ESCIMO_STEP(my_time, height, step_time)
  Symbolics.@register_symbolic ESCIMO_ZIDZ(A, B)
  Symbolics.@register_symbolic ESCIMO_ln(x)
  Symbolics.@register_symbolic Modelica_Blocks_Tables_Internal_getDerTable1DValueNoDer(tableID::Ptr{Nothing}, icol, u, der_u)
  Symbolics.@register_symbolic Modelica_Blocks_Tables_Internal_getTable1DAbscissaUmax(tableID::Ptr{Nothing})
  Symbolics.@register_symbolic Modelica_Blocks_Tables_Internal_getTable1DAbscissaUmin(tableID::Ptr{Nothing})
  Symbolics.@register_symbolic Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(tableID::Ptr{Nothing}, icol, u)
  Symbolics.@register_symbolic Modelica_Blocks_Types_ExternalCombiTable1D_constructor(tableName, fileName, table, columns, smoothness, extrapolation, verboseRead)
  Symbolics.@register_symbolic Modelica_Blocks_Types_ExternalCombiTable1D_destructor(externalCombiTable1D::Ptr{Nothing})
  Symbolics.@register_symbolic Modelica_Utilities_Strings_Advanced_skipWhiteSpace(string, startIndex)
  Symbolics.@register_symbolic Modelica_Utilities_Strings_isEmpty(string)
  Symbolics.@register_symbolic Modelica_Utilities_Strings_length(string)
  begin
    begin
      saved_values_ESCIMO = SavedValues(Float64, Tuple{Float64,Array})
      function ESCIMOCallbackSet(aux)
        local p = aux[1]
        local reals = aux[2]
        local reducedSystem = aux[3]
        begin
          condition4 = ((x, t, integrator) -> t - 2.0e7)
          affect4! = (integrator -> begin

            local t = integrator.t + integrator.dt
            local x = integrator.u
            if integrator.dt == 0.0
              println("integrator.dt was zero. Aborting.")
            end

            if Bool(t >= 2.0e7)
              x[118] = x[790]
            end
          end)
          cb4 = ContinuousCallback(condition4, affect4!, rootfind = true, save_positions = (true, true), affect_neg! = affect4!)
          nothing
        end
        begin
          condition6 = ((x, t, integrator) -> t - 3.0e6)
          affect6! = (integrator -> begin

            local t = integrator.t + integrator.dt
            local x = integrator.u
            if integrator.dt == 0.0
              println("integrator.dt was zero. Aborting.")
            end

            if Bool(t >= 3.0e6)
              x[120] = x[52] * 0.75
            end
          end)
          cb6 = ContinuousCallback(condition6, affect6!, rootfind = true, save_positions = (true, true), affect_neg! = affect6!)
          nothing
        end
        begin
          condition7 = ((x, t, integrator) -> t - 1970.0)
          affect7! = (integrator -> begin

            local t = integrator.t + integrator.dt
            local x = integrator.u
            if integrator.dt == 0.0
              println("integrator.dt was zero. Aborting.")
            end

            if Bool(t >= 1970.0)
                      x[29] = x[741]
                      x[31] = x[795]
            end
          end)
          cb7 = ContinuousCallback(condition7, affect7!, rootfind = true, save_positions = (true, true), affect_neg! = affect7!)
          nothing
        end
        begin
          condition10 = ((x, t, integrator) -> t - 2010.0)
          affect10! = (integrator -> begin

            local t = integrator.t + integrator.dt
            local x = integrator.u
            if integrator.dt == 0.0
              println("integrator.dt was zero. Aborting.")
            end

                    if Bool(t >= 2010.0)
                       x[24] = x[411]
                       x[25] = x[617]
                       x[32] = x[801]
                       x[33] = x[802]
                       x[30] = x[748]
                       x[27] = x[604]
            end
          end)
          cb10 = ContinuousCallback(condition10, affect10!, rootfind = true, save_positions = (true, true), affect_neg! = affect10!)
          nothing
        end
        begin
          condition12 = ((x, t, integrator) -> t - 500000.0)
          affect12! = (integrator -> begin

            local t = integrator.t + integrator.dt
            local x = integrator.u
            if integrator.dt == 0.0
              println("integrator.dt was zero. Aborting.")
            end
            if Bool(t >= 500000.0)
              x[126] = x[737]
            end
          end)
          cb12 = ContinuousCallback(condition12, affect12!, rootfind = true, save_positions = (true, true), affect_neg! = affect12!)
          nothing
        end
        nothing
        return CallbackSet(cb4, cb6, cb7, cb10, cb12)
      end
    end
    function ESCIMOModel(tspan = (0.0, 1.0))
      ModelingToolkit.@variables t
      D = ModelingToolkit.Differential(t)
      parameters = ModelingToolkit.@parameters(begin
        Future_volcanic_emissions
        Albedo_Antarctic_hist
        Albedo_Antarctic_sens
        Albedo_BARREN_normal
        Albedo_BARREN_white
        Albedo_DESERT_normal
        Albedo_glacier_hist
        Albedo_glacier_sens
        Albedo_GRASS_burnt
        Albedo_GRASS_deforested
        Albedo_GRASS_normal_cover
        Albedo_Greenland
        Albedo_NF_burnt
        Albedo_NF_deforested
        Albedo_NF_normal_cover
        Albedo_TROP_burnt
        Albedo_TROP_deforested
        Albedo_TROP_normal_cover
        Albedo_TUNDRA_burnt
        Albedo_TUNDRA_deforested
        Albedo_TUNDRA_normal_cover
        Albedo_URBAN_normal
        Amount_methane_hydrates__clathrates__experimentally_released_GtC
        Amt_of_constant_emissions_GtC_yr
        Annual_pct_increase_CH4_emissions_from_2015_pct_yr
        Annual_pct_increase_CO2_emissions_from_2015_pct_yr
        Antarctic_ice_volume_in_1850_km3
        Arctic_ice_albedo_1850
        Arctic_ice_area_in_1850_km2
        Arctic_surface_temp_delay_yr
        Area_covered_by_high_clouds_in_1850
        Area_covered_by_low_clouds_in_1850
        Area_equivalent_of_1km_linear_retreat_km2
        Area_of_earth_m2
        Area_of_ocean_at_surface_361900_Gm2
        Atmos_heat_used_for_melting_Initially_1_yr
        Average_thickness_arctic_ice_km
        Avg_amount_of_C_in_the_form_of_CH4_per_km2
        Avg_depth_of_permafrost_km
        Avg_flatness_of_worlds_coastline
        Avg_thickness_Antarctic_hist_km
        Avg_thickness_Antarctic_sens_km
        Avg_thickness_Greenland_km
        C_in_atmosphere_in_1850_GtC
        C_in_the_form_of_CH4_in_atm_1850
        Carbon_per_biomass_tC_per_tBiomass
        CC_in_cold_ocean_0_to_100m_1850_ymoles_per_litre
        CC_in_cold_ocean_downwelling_100m_bottom_1850_ymoles_per_litre
        CC_in_ocean_upwelling_100m_to_1km_1850_ymoles_per_litre
        CC_in_warm_ocean_0_to_100m_1850_ymoles_per_litre
        CC_ocean_deep_1km_to_bottom_1850_ymoles_per_litre
        CH4_concentration_in_2010_ppb
        CH4_halflife_in_atmosphere
        Cold_dense_water_sinking_in_Sverdrup_in_1850
        Constant_anthropogenic_CH4_emissions
        Convection_as_f_of_incoming_solar_in_1850
        conversion_factor_CH4_Gt_to_ppb
        Conversion_from_Kyoto_Flour_amount_to_concentration_ppt_kt
        Conversion_from_Montreal_gases_amount_to_concentration_ppt_kt
        Conversion_Millionkm2_to_km2_Mkm2_km2
        Conversion_of_anthro_aerosol_emissions_to_forcing
        Conversion_of_volcanic_aerosol_emissions_to_CO2_emissions_GtC_pr_VAE
        Conversion_of_volcanic_aerosol_forcing_to_volcanic_aerosol_emissions
        Conversion_ymoles_per_kg_to_pCO2_yatm
        Densitiy_of_water_relative_to_ice
        Duration_of_destruction_yr
        Emissions_of_natural_CH4_GtC_yr
        Emissivity_atm
        Emissivity_surface
        Evaporation_as_fraction_of_incoming_solar_in_1850
        EXP_12f_Stratospheric_scattering_experiment_0_off_1_on
        Experimental_doubling_of_constant_C_emissions_how_long_yr
        Experimental_doubling_of_constant_C_emissions_how_much_1_100pct
        Experimental_doubling_of_constant_C_emissions_when_yr
        Frac_of_surface_emission_through_atm_window
        Frac_SW_clear_sky_reflection_aka_scattering
        Frac_SW_HI_cloud_efffect_aka_cloud_albedo
        Frac_SW_LO_cloud_efffect_aka_cloud_albedo
        Fraction_of_C_released_from_permafrost_released_as_CH4_hist_dmnl
        Fraction_of_C_released_from_permafrost_released_as_CH4_sensitivity_dmnl
        Fraction_of_earth_surface_as_ocean
        Fraction_of_heat_needed_to_melt_antarctic_ice_coming_from_air
        Fraction_of_heat_needed_to_melt_arctic_ice_coming_from_air
        Fraction_of_heat_needed_to_melt_Greenland_ice_that_slid_into_the_ocean_coming_from_air
        Fraction_of_methane_hydrates_released_from_the_ocean_converted_to_CO2_before_it_is_relased_to_the_atmosphere
        Fraction_of_ocean_classified_warm_surface
        Glacial_ice_volume_in_1850_km3
        Global_Warming_Potential_CH4
        Global_Warming_Potential_N20
        GRASS_area_burned_in_1850_Mkm2
        GRASS_area_deforested_in_1850_Mkm2
        GRASS_area_harvested_in_1850_Mkm2
        GRASS_Avg_life_biomass_yr
        GRASS_Avg_life_of_building_yr
        GRASS_Biomass_locked_in_construction_material_in_1850_GtBiomass
        GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_in_1850_GtBiomass
        GRASS_Fraction_of_construction_waste_burned_0_1
        GRASS_fraction_of_DeadB_and_SOM_being_destroyed_by_deforesting
        GRASS_fraction_of_DeadB_and_SOM_being_destroyed_by_energy_harvesting
        GRASS_fraction_of_DeadB_and_SOM_destroyed_by_natural_fires
        GRASS_living_biomass_densitiy_in_1850_tBiomass_pr_km2
        GRASS_Living_biomass_in_1850_GtBiomass
        GRASS_Normal_fire_incidence_fraction_yr
        GRASS_Ref_historical_deforestation_pct_yr
        GRASS_runoff_time
        GRASS_Speed_of_regrowth_yr
        GRASS_Time_to_decompose_undisturbed_dead_biomass_yr
        Greenland_ice_slide_circulation_slowdown_effect
        Greenland_ice_volume_in_1850_km3
        Greenland_slide_experiment_how_much_sildes_in_the_ocean_fraction
        Greenland_slide_experiment_over_how_many_years_yr
        GtIce_vs_km3
        Heat_gained___needed_to_freeze___unfreeze_1_km3_permafrost_ZJ_km3
        Heat_in__ocean__deep_in_1850_ZJ
        Heat_in_atmosphere_in_1850_ZJ
        Heat_in_surface_in_1850_ZJ
        Heat_needed_to_melt_1_km3_of_ice_ZJ
        Hist_Avg_thickness_glacier_km
        Hist_Net_heat_flow_ocean_between_surface_and_deep_per_K_of_difference_ZJ_yr_K
        Hist_NF_Avg_life_biomass_yr
        Hist_NF_Speed_of_regrowth_yr
        Hist_Slope_of_effect_of_temp_shifting_DESERT_to_GRASS
        Hist_Slope_temp_vs_glacial_ice_melting
        Hist_Time_in_trunk
        Hist_Time_to_degrade_Kyoto_Flour_yr
        Hist_Time_to_regrow_NF_after_buning_yr
        Hist_TROP_runoff_time
        Hist_TROP_Time_to_decompose_undisturbed_dead_biomass_yr
        K_to_C_conversion_C_K
        Kyoto_Flour_Global_Warming_Potential
        Land_surface_temp_adjustment_time_yr
        LW_ALL_cloud_radiation_reference_in_1850_W_m2
        LW_LO_cloud_radiation_reference_in_1850_W_m2
        LW_radiation_fraction_blocked_by_other_GHG_in_1850
        Man_made_CH4_emissions_in_2015_GtC
        Man_made_CO2_emissions_in_2015_GtC
        MAX_NATURE_CCS_removal_in_2050_GtCO2e_yr
        Melting_of_permafrost_at_all_depths_at_4_deg_C_temp_diff_km_yr
        Montreal_Global_Warming_Potential
        Myhre_constant_for_CH4
        Myhre_constant_for_CO2
        Myhre_constant_for_N20
        N2O_concentration_in_2010_ppb
        N2O_in_atmosphere_MtN2O_in_1850
        N2O_natural_emissions
        Net_marine_primary_production_in_1850
        NEvt_13a_double_rate_of_melting_ice_and_permafrost
        NEvt_13b2_Double_incidence_of_biomass_fires
        NEvt_13b3_double_sunspot_amplitude_from_2015_onwards_1_normal_2_double
        NEvt_13c1_increase_in_area_covered_by_low_clouds
        NEvt_13d_Greenland_slide_experiment_start_yr
        NEvt_2a_Volcanic_eruptions_in_the_future_VAEs_first_future_pulse
        NEvt_3b_increase_in_area_covered_by_high_clouds
        NF_area_burned_in_1850_Mkm2
        NF_area_deforested_in_1850_Mkm2
        NF_area_harvested_in_1850_Mkm2
        NF_Avg_life_of_building_yr
        NF_Biomass_locked_in_construction_material_in_1850_GtBiomass
        NF_Dead_biomass__litter_and_soil_organic_matter_SOM_in_1850_GtBiomass
        NF_DeadB_and_SOM_densitiy_in_1850_tBiomass_pr_km2
        NF_Fraction_of_construction_waste_burned_0_1
        NF_fraction_of_DeadB_and_SOM_being_destroyed_by_clear_cutting
        NF_fraction_of_DeadB_and_SOM_being_destroyed_by_deforesting
        NF_fraction_of_DeadB_and_SOM_being_destroyed_by_energy_harvesting
        NF_fraction_of_DeadB_and_SOM_destroyed_by_natural_fires
        NF_living_biomass_densitiy_in_1850_tBiomass_pr_km2
        NF_Living_biomass_in_1850_GtBiomass
        NF_Normal_fire_incidence_fraction_yr
        NF_Ref_historical_deforestation_pct_yr
        NF_runoff_time
        NF_Time_to_decompose_undisturbed_dead_biomass_yr
        Ocean_heat_used_for_melting_Initially_1_yr
        Ocean_slowdown_experimental_factor
        Open_ocean_albedo
        Over_how_many_yrs_methane_hydrate_release_yr
        per_annum_yr
        Policy_1_Reducing_GHG_emissions_by_one_third_by_2035
        Policy_2_Large_scale_implementation_of_carbon_capture_and_geological_storage__CCS_
        Population_2000_bn
        Pressure_adjustment_deep_pct
        Pressure_adjustment_surface_pct
        Rate_of_wetland_destruction_pct_of_existing_wetlands_yr
        Ratio_of_methane_in_tundra_to_wetland
        Ref_shifting_biome_yr
        Ref_temp_difference__4_degC_
        Ref_temp_difference_for_antarctic_ice_melting__3_degC_
        Ref_temp_difference_for_Arctic_ice_melting
        Ref_temp_difference_for_glacial_ice_melting__1_degC_
        Ref_temp_difference_for_greenland_ice_melting_C
        Ref_temp_difference_for_greenland_ice_that_slid_into_the_ocean_melting_degC
        Reference_temp_C
        Reference_Time_to_regrow_TROP_after_deforesting_yr
        SCALE_and_UNIT_converter_zero_C_to_K
        Sens_Avg_thickness_glacier_km
        Sens_Frac_atm_absorption
        Sens_Net_heat_flow_ocean_between_surface_and_deep_per_K_of_difference_ZJ_yr_K
        Sens_NF_Avg_life_biomass_yr
        Sens_NF_Speed_of_regrowth_yr
        Sens_Slope_of_effect_of_temp_shifting_DESERT_to_GRASS
        Sens_Slope_temp_vs_glacial_ice_melting
        Sens_Time_in_trunk
        Sens_Time_to_degrade_Kyoto_Flour_yr
        Sens_Time_to_regrow_NF_after_buning_yr
        Sens_TROP_runoff_time
        Sens_TROP_Time_to_decompose_undisturbed_dead_biomass_yr
        Sensitivity_of_biomass_new_growth_to_CO2_concentration
        Sensitivity_of_convection_to_temp
        Sensitivity_of_evaporation_to_temp
        Sensitivity_of_high_cloud_coverage_to_temp_base
        Sensitivity_of_high_cloud_coverage_to_temp_sens
        Sensitivity_of_low_cloud_coverage_to_temp
        Sensitivity_of_trop_to_humidity
        Slider_for_annual_removal_of_C_from_atm_after_2020_GtC_y
        Slider_for_H2O_slope_hist
        Slider_for_slope_fut
        Slope_btw_Kyoto_Flour_ppt_and_blocking_multiplier
        Slope_btw_Montreal_gases_ppt_and_blocking_multiplier
        Slope_btw_N2O_ppb_and_blocking_multiplier
        Slope_btw_temp_and_permafrost_melting___freezing_base
        Slope_btw_temp_and_permafrost_melting___freezing_sensitivity
        Slope_Effect_Temp_on_NMPP
        Slope_of_effect_of_temp_on_shifting_NF_to_Tundra
        Slope_of_effect_of_temp_on_shifting_TROP_to_NF
        Slope_of_effect_of_temp_shifting_GRASS_to_DESERT
        Slope_of_effect_of_temp_shifting_GRASS_to_NF
        Slope_of_effect_of_temp_shifting_GRASS_to_TROP
        Slope_of_effect_of_temp_shifting_NF_to_GRASS
        Slope_of_effect_of_temp_shifting_NF_to_TROP
        Slope_of_effect_of_temp_shifting_TROP_to_GRASS
        Slope_of_effect_of_temp_shifting_tundra_to_NF
        Slope_of_efffect_of_acidification_on_NMPP
        Slope_temp_eff_on_fire_incidence
        Slope_temp_vs_antarctic_ice_melting
        Slope_temp_vs_Arctic_ice_melting
        Slope_temp_vs_greenland_ice_melting
        Slope_temp_vs_greenland_ice_that_slid_into_the_ocean_melting
        Solar_sine_forcing_amplitude
        Solar_sine_forcing_lift
        Solar_sine_forcing_offset_yr
        Solar_sine_forcing_period_yr
        Stephan_Boltzmann_constant
        Stratospheric_scattering_experiment_end_year
        Stratospheric_scattering_experiment_reduction_from_2015_in_W_m2
        Switch_0_normal_model_1_dbl_CO2_2_1pct_incr
        Switch_btw_historical_CO2_CH4_emissions_or_constant_1history_0constant
        SWITCH_for_NATURE_comm_200115_base_1_cut_all_mm_emi_in_2020_2
        SWITCH_future_slope_base_0_plus_5_1_minus_5_2
        SWITCH_h2o_blocked_table_0_linear_1_poly_2
        SWITCH_h2o_poly_dyn_0_equ_1
        SWITCH_nature_rev_0_base_1_steeper_2_less_steep
        Switch_to_choose_CO2_CH4_aerosol_other_GHG_0normal_1zero_from_2010_2constant_from_2010
        Switch_to_choose_input_emission_scenario_for_CO2_CH4_and_oth_GHG
        Switch_to_drive_model_with_normal_ESCIMO_data__0__CO2e_from_C_Roads__1__or_CO2e_from_CAT_2__or_user_determined_CO2_max_to_find_temp_tipping_point__3_
        Switch_to_run_experiment_12a_reduction_in_emissions_0_off_1_on
        Switch_to_run_experiment_12b_CCS_0_off_1_on
        Switch_to_run_experiment_12c_stopping_TROP_deforestation_0_off_1_on
        Switch_to_run_experiment_12e_white_surfaces_0_off_1_on
        Switch_to_run_NATURE_experiment_CCS_0_off_1_on_0
        Switch_to_run_POLICY_4_Stopping_logging_in_Northern_forests_0_off_1_on
        Temp__ocean__deep_in_1850_C
        Temp_atm_1850
        Temp_gradient_in_surface_degK
        Temp_surface_1850_K
        TEST_Year_in_which_zero_emissions_are_to_be_reached_yr_Remember_to_set_switch_to_9Linear
        Thickness_of_deep_water_box_1km_to_bottom
        Thickness_of_intermediate_water_box_800m
        Thickness_of_surface_water_box_100m
        Time_at_which_human_deforestation_is_stopped
        Time_for_volcanic_aerosols_to_remain_in_the_stratosphere
        Time_in_cold
        Time_in_deep
        Time_in_intermediate_yr
        Time_in_warm
        Time_to_degrade_Montreal_gases_yr
        Time_to_degrade_N2O_in_atmopshere_yr
        Time_to_deposit_C_in_sediment
        Time_to_let_shells_form_and_sink_to_sediment_yr
        Time_to_melt_Arctic_ice_at_the_reference_delta_temp
        Time_to_melt_greenland_ice_at_the_reference_delta_temp
        Time_to_melt_greenland_ice_that_slid_into_the_ocean_at_the_reference_delta_temp
        Time_to_melt_or_freeze_antarctic_ice_at_the_reference_delta_temp
        Time_to_melt_or_freeze_glacial_ice_at_the_reference_delta_temp
        Time_to_propagate_temperature_change_through_the_volume_of_permafrost_yr
        Time_to_reach_C_equilibrium_between_atmosphere_and_ocean
        Time_to_regrow_GRASS_after_buning_yr
        Time_to_regrow_GRASS_after_deforesting_yr
        Time_to_regrow_NF_after_deforesting_yr
        Time_to_regrow_TROP_after_buning_yr
        Time_to_regrow_TUNDRA_after_buning_yr
        Time_to_regrow_TUNDRA_after_deforesting_yr
        Time_to_smooth_out_temperature_diff_relevant_for_melting_or_freezing_from_1850_yr
        Tipping_point_search_amount_at_peak
        Tipping_point_year_of_end
        Tipping_point_year_of_start
        TROP_area_burned_in_1850_Mkm2
        TROP_area_deforested_in_1850_Mkm2
        TROP_area_harvested_in_1850_Mkm2
        TROP_Avg_life_biomass_yr
        TROP_Avg_life_of_building_yr
        TROP_Biomass_locked_in_construction_material_in_1850_GtBiomass
        TROP_clear_cut_fraction
        TROP_Dead_biomass__litter_and_soil_organic_matter_SOM_in_1850_GtBiomass
        TROP_DeadB_and_SOM_densitiy_in_1850_tBiomass_pr_km2
        TROP_Fraction_of_construction_waste_burned_0_1
        TROP_fraction_of_DeadB_and_SOM_being_destroyed_by_clear_cutting
        TROP_fraction_of_DeadB_and_SOM_being_destroyed_by_deforesting
        TROP_fraction_of_DeadB_and_SOM_being_destroyed_by_energy_harvesting
        TROP_fraction_of_DeadB_and_SOM_destroyed_by_natural_fires
        TROP_living_biomass_densitiy_in_1850_tBiomass_pr_km2
        TROP_Living_biomass_in_1850_GtBiomass
        TROP_Normal_fire_incidence_fraction_yr
        TROP_Ref_historical_deforestation_pct_yr
        TROP_Slope_temp_eff_on_potential_biomass_per_km2
        TROP_Speed_of_regrowth_yr
        TUNDRA_area_burned_in_1850_Mkm2
        TUNDRA_area_deforested_in_1850_Mkm2
        TUNDRA_area_harvested_in_1850_Mkm2
        TUNDRA_Avg_life_biomass_yr
        TUNDRA_Avg_life_of_building_yr
        TUNDRA_Biomass_locked_in_construction_material_in_1850_GtBiomass
        TUNDRA_Dead_biomass__litter_and_soil_organic_matter_SOM_in_1850_GtBiomass
        TUNDRA_DeadB_and_SOM_densitiy_in_1850_tBiomass_pr_km2
        TUNDRA_Fraction_of_construction_waste_burned_0_1
        TUNDRA_fraction_of_DeadB_and_SOM_being_destroyed_by_deforesting
        TUNDRA_fraction_of_DeadB_and_SOM_being_destroyed_by_energy_harvesting
        TUNDRA_fraction_of_DeadB_and_SOM_destroyed_by_natural_fires
        TUNDRA_living_biomass_densitiy_in_1850_tBiomass_pr_km2
        TUNDRA_Living_biomass_in_1850_GtBiomass
        TUNDRA_Normal_fire_incidence_fraction_yr
        TUNDRA_Ref_historical_deforestation_pct_yr
        TUNDRA_runoff_time
        TUNDRA_Speed_of_regrowth_yr
        TUNDRA_Time_to_decompose_undisturbed_dead_biomass_yr
        UNIT_conversion_1_km3
        UNIT_conversion_1_yr
        UNIT_conversion_C_to_pH
        UNIT_Conversion_from__km3__km_yr___to_Mkm2_yr
        UNIT_conversion_from_km_to_m
        UNIT_Conversion_from_km3_to_km2
        UNIT_Conversion_from_N2O_amount_to_concentration_ppb_MtN2O
        UNIT_conversion_Gm3_to_km3
        UNIT_conversion_Gt_to_kt
        UNIT_conversion_Gt_to_Mt
        UNIT_conversion_GtBiomass_yr_to_Mkm2_yr
        UNIT_conversion_GtC_to_MtC
        UNIT_conversion_GtIce_to_ZJ_melting
        UNIT_conversion_km2___km_to_km3
        UNIT_conversion_km2_to_Mkm2
        UNIT_conversion_km3_to_Gm3
        UNIT_conversion_km3_km_to_km2
        UNIT_conversion_m2_to_km2
        UNIT_conversion_m2_to_Mkm2
        UNIT_conversion_Sv_to_Gm3_yr
        UNIT_conversion_to_Gm3
        UNIT_conversion_to_km2_yr
        UNIT_conversion_to_yr
        UNIT_conversion_W_to_ZJ_s
        UNIT_conversion_ymoles___litre_to_dless
        UNIT_conversion_yr_to_dless
        Urban_area_fraction_2000
        Use_of_GRASS_biomass_for_construction_in_1850_pct
        Use_of_GRASS_biomass_for_energy_in_1850_pct
        Use_of_NF_biomass_for_construction_in_1850_pct
        Use_of_NF_biomass_for_energy_in_1850_pct
        Use_of_TROP_biomass_for_construction_in_1850_pct
        Use_of_TROP_biomass_for_energy_in_1850_pct
        Use_of_TUNDRA_biomass_for_construction_in_1850_pct
        Use_of_TUNDRA_biomass_for_energy_in_1850_pct
        VAES_puls_repetition
        VAES_pulse_duration
        VAES_pulse_height
        Value_of_anthropogenic_aerosol_emissions_during_2015
        Water_content_of_evaporation_g_kg_per_ZJ_yr
        Wetlands_area_1850
        When_first_destroyed_yr
        When_methane_hydrates_first_released_yr
        When_to_sample_for_CO2_experiment_yr
        Yr_to_cut_mm_emi_abrubtly_to_zero_y
        Zero_C_on_K_scale_K
        Zetta
        CO2_concentration_in_1750_ppm
        N2O_ie_N_1750_ppb
        CH4_ie_M_1750_ppb
        LW_Clear_sky_emissions_from_atm_W_m2_in_1850
        SW_surface_absorption_W_m2_in_1850
        SW_surface_reflection_W_m2_in_1850
        C_in_TUNDRA_DeadB_and_soil_in_1850_GtC
        C_in_TUNDRA_LB_in_1850_GtC
        Ga__BB_radiation_less_TOA_radiation_W_m2_in_1850
        Biomass_new_growing_1850_GtBiomass___yr
        Switch_to_choose_CO2_CH4_aerosol_other_GHG_0normal_1zero_from_202constant_from_2010
        LW_TOA_radiation_from_atm_to_space_in_1850_W_m2
      end)
      begin
        variableConstructors = Function[]
        begin
          function generateStateVariables1()
            (
              :ifCond1,
              :ifCond2,
              :Model_N2O_concentration_in_1850_ppb,
              :CO2_concentration_in_1850_ppm,
              :Incoming_solar_in_1850_ZJ_yr,
              :C_in_atmosphere_GtC_in_1850,
              :C_in_biomass_in_1850_GtC,
              :Total_carbon_in_ocean_GtC_in_1850,
              :Temp_ocean_deep_1850_degC,
              :init_ph_in_cold_water,
              :Humidity_of_atmosphere_in_1850_g_kg,
              :LW_TOA_radiation_from_atm_to_space_in_1850,
              :Temp__ocean__surface_in_1850_C,
              :Fraction_blocked_by_ALL_GHG_in_1850,
              :Fraction_blocked_CO2_in_1850,
              :Fraction_blocked_CH4_in_1850,
              :Fraction_blocked_othGHG_in_1850,
              :init_C_in_GRASS,
              :init_C_in_NF,
              :init_C_in_TROP,
              :init_C_in_TUNDRA,
              :Fossil_fuel_reserves_in_ground_1850_GtC,
              :Time,
              :Aerosol_anthropogenic_emissions_in_2010,
              :CO2_emissions_in_2010,
              :CO2_ppm_value_at_When_to_sample,
              :CO4_emissions_in_2010,
              :Greenland_slide_experiment_end_condition,
              :Kyoto_Flour_concentration_in_1970_ppt,
              :Kyoto_Flour_emissions_RCPs_JR_in_2010,
              :Montreal_gases_concentration_in_1970_ppt,
              :Montreal_gases_emissions_RCPs_JR_in_2010,
              :N20_emissions_RCPs_JR_in_2010,
              :Tipping_point_search_amount_at_start,
              :Arctic_land_surface_temp_anomaly_compared_to_1850,
              :Biological_removal_of_C_from_WSW_GtC_per_yr,
              :Effect_of_temp_on_permafrost_melting_dmnl,
              :Temp_diff_relevant_for_melting_or_freezing_arctic_ice_from_1850,
              :Temp_diff_relevant_for_melting_or_freezing_from_1850,
              :yr_on_yr_change_in_C_in_atm_GtC_yr,
              :C_in_ocean_1_yr_ago_GtC,
              :C_in_ocean_1_yr_ago_GtC_LV1,
              :C_in_ocean_1_yr_ago_GtC_LV2,
              :Atmos_heat_used_for_melting_last_year_1_yr_LV,
              :Ocean_heat_used_for_melting_last_year_ZJ_yr_LV,
              :C_in_atm_1_yr_ago_GtC_LV1,
              :C_in_atm_1_yr_ago_GtC_LV2,
              :C_in_atm_1_yr_ago_GtC_LV3,
              :All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_LV1,
              :All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_LV2,
            )
          end
          push!(variableConstructors, generateStateVariables1)
        end
        begin
          function generateStateVariables2()
            (
              :All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_LV3,
              :Antarctic_ice_volume_km3,
              :Arctic_ice__on_sea__area_km2,
              :C_in_atmosphere_GtC,
              :C_in_atmosphere_in_form_of_CH4,
              :C_in_cold_surface_water_GtC,
              :C_in_cold_water_trunk_downwelling_GtC,
              :C_in_deep_water_volume_1km_to_bottom_GtC,
              :C_in_intermediate_upwelling_water_100m_to_1km_GtC,
              :C_in_permafrost_in_form_of_CH4,
              :C_in_sediment,
              :C_in_warm_surface_water_GtC,
              :Cold_surface_water_volume_Gm3,
              :Cold_water_volume_downwelling_Gm3,
              :Cumulative_antarctic_ice_volume_loss_GtIce,
              :Cumulative_C_released_from_permafrost_as_either_CH4_or_CO2_in_GtC_yr,
              :Cumulative_carbon_captured_and_stored_GtC,
              :Cumulative_carbon_removed_from_atm_for_nature_May_2020,
              :Cumulative_flow_of_C_to_biomass_since_1850_GtC,
              :Cumulative_glacial_ice_volume_loss_GtIce,
              :Cumulative_Greenland_ice_volume_loss_GtIce,
              :Cumulative_heat_to_atm_ZJ,
              :Cumulative_ocean_volume_increase_due_to_ice_melting_km3,
              :Cumulative_release_of_C_from_permafrost_GtC,
              :Deep_water_volume_1km_to_4km_Gm3,
              :DESERT_Mkm2,
              :Fossil_fuel_reserves_in_ground_GtC,
              :Glacial_ice_volume_km3,
              :GRASS_area_burnt_Mkm2,
              :GRASS_area_harvested_Mkm2,
              :GRASS_Biomass_locked_in_construction_material_GtBiomass,
              :GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass,
              :GRASS_deforested_Mkm2,
              :GRASS_Living_biomass_GtBiomass,
              :GRASS_potential_area_Mkm2,
              :Greenland_ice_volume_on_Greenland_km3,
              :Greenland_ice_volume_that_slid_into_the_ocean_km3,
              :Heat_in_atmosphere_ZJ,
              :Heat_in_deep_ZJ,
              :Heat_in_surface,
              :Intermediate_upwelling_water_volume_100m_to_1km_Gm3,
              :Kyoto_Flour_gases_in_atm,
              :Montreal_gases_in_atm,
              :N2O_in_atmosphere_MtN2O,
              :NATURE_Cumulative_CCS_GtC,
              :NF_area_burnt_Mkm2,
              :NF_area_clear_cut_Mkm2,
              :NF_area_deforested_Mkm2,
              :NF_area_harvested_Mkm2,
              :NF_Biomass_locked_in_construction_material_GtBiomass,
            )
          end
          push!(variableConstructors, generateStateVariables2)
        end
        begin
          function generateStateVariables3()
            (
              :NF_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass,
              :NF_Living_biomass_GtBiomass,
              :NF_potential_area_Mkm2,
              :Sum_C_absorbed_by_ocean_GtC,
              :Sum_heat_to_deep_ocean,
              :Sum_heat_to_deep_ocean_btw_72_and_08,
              :Sum_heat_to_surface_ocean_btw_72_and_08,
              :Sum_heat_to_surface_ocean_ZJ,
              :Sum_man_made_CO2_emissions_GtC,
              :Sum_net_C_to_atm,
              :TROP_area_burnt_Mkm2,
              :TROP_area_clear_cut_Mkm2,
              :TROP_area_deforested_Mkm2,
              :TROP_area_harvested_Mkm2,
              :TROP_Biomass_locked_in_construction_material_GtBiomass,
              :TROP_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass,
              :TROP_Living_biomass_GtBiomass,
              :TROP_potential_area_Mkm2,
              :TUNDRA_area_burnt_Mkm2,
              :TUNDRA_area_harvested_Mkm2,
              :TUNDRA_Biomass_locked_in_construction_material_GtBiomass,
              :TUNDRA_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass,
              :TUNDRA_deforested_Mkm2,
              :TUNDRA_Living_biomass_GtBiomass,
              :Tundra_potential_area_Mkm2,
              :Volcanic_aerosols_in_stratosphere,
              :Warm_surface_water_volume_Gm3,
              :Wetlands_area,
            )
          end
          push!(variableConstructors, generateStateVariables3)
        end
        begin
          function generateAlgebraicVariables1()
            (
              :combi_E3_SC_1_CO2_GtC_yr_u,
              Symbol("combi_E3_SC_1_CO2_GtC_yr_y[1]"),
              :E3_SC_1_CO2_GtC_yr,
              :combi_E3_SC_1_CH4_GtC_yr_u,
              Symbol("combi_E3_SC_1_CH4_GtC_yr_y[1]"),
              :E3_SC_1_CH4_GtC_yr,
              :combi_E3_SC_1_N2O_Mt_yr_u,
              Symbol("combi_E3_SC_1_N2O_Mt_yr_y[1]"),
              :E3_SC_1_N2O_Mt_yr,
              :combi_E3_SC_1_Kyoto_F_kt_yr_u,
              Symbol("combi_E3_SC_1_Kyoto_F_kt_yr_y[1]"),
              :E3_SC_1_Kyoto_F_kt_yr,
              :combi_E3_SC_1_Montreal_gases_kt_yr_u,
              Symbol("combi_E3_SC_1_Montreal_gases_kt_yr_y[1]"),
              :E3_SC_1_Montreal_gases_kt_yr,
              :combi_E3_SC_2_CO2_GtC_yr_u,
              Symbol("combi_E3_SC_2_CO2_GtC_yr_y[1]"),
              :E3_SC_2_CO2_GtC_yr,
              :combi_E3_SC_2_CH4_GtC_yr_u,
              Symbol("combi_E3_SC_2_CH4_GtC_yr_y[1]"),
              :E3_SC_2_CH4_GtC_yr,
              :combi_E3_SC_2_N2O_Mt_yr_u,
              Symbol("combi_E3_SC_2_N2O_Mt_yr_y[1]"),
              :E3_SC_2_N2O_Mt_yr,
              :combi_E3_SC_2_Kyoto_F_kt_yr_u,
              Symbol("combi_E3_SC_2_Kyoto_F_kt_yr_y[1]"),
              :E3_SC_2_Kyoto_F_kt_yr,
              :combi_E3_SC_2_Montreal_gases_kt_yr_u,
              Symbol("combi_E3_SC_2_Montreal_gases_kt_yr_y[1]"),
              :E3_SC_2_Montreal_gases_kt_yr,
              :combi_E3_SC_3_CO2_GtC_yr_u,
              Symbol("combi_E3_SC_3_CO2_GtC_yr_y[1]"),
              :E3_SC_3_CO2_GtC_yr,
              :combi_E3_SC_3_CH4_GtC_yr_u,
              Symbol("combi_E3_SC_3_CH4_GtC_yr_y[1]"),
              :E3_SC_3_CH4_GtC_yr,
              :combi_E3_SC_3_N2O_Mt_yr_u,
              Symbol("combi_E3_SC_3_N2O_Mt_yr_y[1]"),
              :E3_SC_3_N2O_Mt_yr,
              :combi_E3_SC_3_Kyoto_F_kt_yr_u,
              Symbol("combi_E3_SC_3_Kyoto_F_kt_yr_y[1]"),
              :E3_SC_3_Kyoto_F_kt_yr,
              :combi_E3_SC_3_Montreal_gases_kt_yr_u,
              Symbol("combi_E3_SC_3_Montreal_gases_kt_yr_y[1]"),
              :E3_SC_3_Montreal_gases_kt_yr,
              :combi_E3_SC_4_CO2_GtC_yr_u,
              Symbol("combi_E3_SC_4_CO2_GtC_yr_y[1]"),
              :E3_SC_4_CO2_GtC_yr,
              :combi_E3_SC_4_CH4_GtC_yr_u,
              Symbol("combi_E3_SC_4_CH4_GtC_yr_y[1]"),
            )
          end
          push!(variableConstructors, generateAlgebraicVariables1)
        end
        begin
          function generateAlgebraicVariables2()
            (
              :E3_SC_4_CH4_GtC_yr,
              :combi_E3_SC_4_N2O_Mt_yr_u,
              Symbol("combi_E3_SC_4_N2O_Mt_yr_y[1]"),
              :E3_SC_4_N2O_Mt_yr,
              :combi_E3_SC_4_Kyoto_F_kt_yr_u,
              Symbol("combi_E3_SC_4_Kyoto_F_kt_yr_y[1]"),
              :E3_SC_4_Kyoto_F_kt_yr,
              :combi_E3_SC_4_Montreal_gases_kt_yr_u,
              Symbol("combi_E3_SC_4_Montreal_gases_kt_yr_y[1]"),
              :E3_SC_4_Montreal_gases_kt_yr,
              :combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_MtCH4_yr_u,
              Symbol("combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_MtCH4_yr_y[1]"),
              :Emissions_of_anthro_CH4_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_MtCH4_yr,
              :combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP3_MtCH4_yr_u,
              Symbol("combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP3_MtCH4_yr_y[1]"),
              :Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP3_MtCH4_yr,
              :combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP45_MtCH4_yr_u,
              Symbol("combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP45_MtCH4_yr_y[1]"),
              :Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP45_MtCH4_yr,
              :combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP6_MtCH4_yr_u,
              Symbol("combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP6_MtCH4_yr_y[1]"),
              :Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP6_MtCH4_yr,
              :combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP85_MtCH4_yr_u,
              Symbol("combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP85_MtCH4_yr_y[1]"),
              :Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP85_MtCH4_yr,
              :combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_GtC_yr_u,
              Symbol("combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_GtC_yr_y[1]"),
              :Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_GtC_yr,
              :combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_GtC_yr_u,
              Symbol("combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_GtC_yr_y[1]"),
              :Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_GtC_yr,
              :combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_hist_with_JR__large_scale_CCS__exp_GtC_yr_u,
              Symbol("combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_hist_with_JR__large_scale_CCS__exp_GtC_yr_y[1]"),
              :Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_hist_with_JR__large_scale_CCS__exp_GtC_yr,
              :combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP3_GtC_yr_u,
              Symbol("combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP3_GtC_yr_y[1]"),
              :Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP3_GtC_yr,
              :combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP45_GtC_yr_u,
              Symbol("combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP45_GtC_yr_y[1]"),
              :Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP45_GtC_yr,
              :combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP6_GtC_yr_u,
              Symbol("combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP6_GtC_yr_y[1]"),
              :Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP6_GtC_yr,
              :combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP85_GtC_yr_u,
              Symbol("combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP85_GtC_yr_y[1]"),
              :Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP85_GtC_yr,
              :combi_Emissions_of_anthro_N20_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_MtN2O_yr_u,
              Symbol("combi_Emissions_of_anthro_N20_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_MtN2O_yr_y[1]"),
              :Emissions_of_anthro_N20_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_MtN2O_yr,
              :combi_Emissions_of_anthro_N20_with_JR_2052_shape_MtN20_yr_u,
            )
          end
          push!(variableConstructors, generateAlgebraicVariables2)
        end
        begin
          function generateAlgebraicVariables3()
            (
              Symbol("combi_Emissions_of_anthro_N20_with_JR_2052_shape_MtN20_yr_y[1]"),
              :Emissions_of_anthro_N20_with_JR_2052_shape_MtN20_yr,
              :combi_Emissions_of_Kyoto_Flour_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_kt_yr_u,
              Symbol("combi_Emissions_of_Kyoto_Flour_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_kt_yr_y[1]"),
              :Emissions_of_Kyoto_Flour_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_kt_yr,
              :combi_Emissions_of_Kyoto_Flour_with_JR_2052_shape_kt_yr_u,
              Symbol("combi_Emissions_of_Kyoto_Flour_with_JR_2052_shape_kt_yr_y[1]"),
              :Emissions_of_Kyoto_Flour_with_JR_2052_shape_kt_yr,
              :combi_Emissions_of_Montreal_gases_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_kt_yr_u,
              Symbol("combi_Emissions_of_Montreal_gases_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_kt_yr_y[1]"),
              :Emissions_of_Montreal_gases_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_kt_yr,
              :combi_Emissions_of_Montreal_gases_with_JR_2052_shape_kt_yr_u,
              Symbol("combi_Emissions_of_Montreal_gases_with_JR_2052_shape_kt_yr_y[1]"),
              :Emissions_of_Montreal_gases_with_JR_2052_shape_kt_yr,
              :combi_CH4_emissions_from_CO2e_C_Roads_u,
              Symbol("combi_CH4_emissions_from_CO2e_C_Roads_y[1]"),
              :CH4_emissions_from_CO2e_C_Roads,
              :combi_CH4_emissions_from_CO2e_CAT_u,
              Symbol("combi_CH4_emissions_from_CO2e_CAT_y[1]"),
              :CH4_emissions_from_CO2e_CAT,
              :combi_CH4_emissions_pct_contribution_to_Total_CO2e_u,
              Symbol("combi_CH4_emissions_pct_contribution_to_Total_CO2e_y[1]"),
              :CH4_emissions_pct_contribution_to_Total_CO2e,
              :combi_CO2_emissions_from_CO2e_C_Roads_u,
              Symbol("combi_CO2_emissions_from_CO2e_C_Roads_y[1]"),
              :CO2_emissions_from_CO2e_C_Roads,
              :combi_CO2_emissions_from_CO2e_CAT_u,
              Symbol("combi_CO2_emissions_from_CO2e_CAT_y[1]"),
              :CO2_emissions_from_CO2e_CAT,
              :combi_CO2_emissions_pct_contribution_to_Total_CO2e_u,
              Symbol("combi_CO2_emissions_pct_contribution_to_Total_CO2e_y[1]"),
              :CO2_emissions_pct_contribution_to_Total_CO2e,
              :combi_Historical_aerosol_emissions_anthro_u,
              Symbol("combi_Historical_aerosol_emissions_anthro_y[1]"),
              :Historical_aerosol_emissions_anthro,
              :combi_Historical_forcing_from_solar_insolation_W_m2_u,
              Symbol("combi_Historical_forcing_from_solar_insolation_W_m2_y[1]"),
              :Historical_forcing_from_solar_insolation_W_m2,
              :combi_Historical_aerosol_forcing_volcanic_u,
              Symbol("combi_Historical_aerosol_forcing_volcanic_y[1]"),
              :Historical_aerosol_forcing_volcanic,
              :combi_OGHG_Kyoto_Flour_emi_rcp3_u,
              Symbol("combi_OGHG_Kyoto_Flour_emi_rcp3_y[1]"),
              :OGHG_Kyoto_Flour_emi_rcp3,
              :combi_OGHG_Kyoto_Flour_emi_rcp45_u,
              Symbol("combi_OGHG_Kyoto_Flour_emi_rcp45_y[1]"),
              :OGHG_Kyoto_Flour_emi_rcp45,
              :combi_OGHG_Kyoto_Flour_emi_rcp6_u,
              Symbol("combi_OGHG_Kyoto_Flour_emi_rcp6_y[1]"),
              :OGHG_Kyoto_Flour_emi_rcp6,
            )
          end
          push!(variableConstructors, generateAlgebraicVariables3)
        end
        begin
          function generateAlgebraicVariables4()
            (
              :combi_OGHG_Kyoto_Flour_emi_rcp85_u,
              Symbol("combi_OGHG_Kyoto_Flour_emi_rcp85_y[1]"),
              :OGHG_Kyoto_Flour_emi_rcp85,
              :combi_Kyoto_Flour_emissions_from_CO2e_C_Roads_u,
              Symbol("combi_Kyoto_Flour_emissions_from_CO2e_C_Roads_y[1]"),
              :Kyoto_Flour_emissions_from_CO2e_C_Roads,
              :combi_Kyoto_Flour_emissions_from_CO2e_CAT_u,
              Symbol("combi_Kyoto_Flour_emissions_from_CO2e_CAT_y[1]"),
              :Kyoto_Flour_emissions_from_CO2e_CAT,
              :combi_Kyoto_Flour_emissions_pct_contribution_to_Total_CO2e_u,
              Symbol("combi_Kyoto_Flour_emissions_pct_contribution_to_Total_CO2e_y[1]"),
              :Kyoto_Flour_emissions_pct_contribution_to_Total_CO2e,
              :combi_OGHG_Montreal_gases_emi_rcp3_u,
              Symbol("combi_OGHG_Montreal_gases_emi_rcp3_y[1]"),
              :OGHG_Montreal_gases_emi_rcp3,
              :combi_OGHG_Montreal_gases_emi_rcp45_u,
              Symbol("combi_OGHG_Montreal_gases_emi_rcp45_y[1]"),
              :OGHG_Montreal_gases_emi_rcp45,
              :combi_OGHG_Montreal_gases_emi_rcp6_u,
              Symbol("combi_OGHG_Montreal_gases_emi_rcp6_y[1]"),
              :OGHG_Montreal_gases_emi_rcp6,
              :combi_OGHG_Montreal_gases_emi_rcp85_u,
              Symbol("combi_OGHG_Montreal_gases_emi_rcp85_y[1]"),
              :OGHG_Montreal_gases_emi_rcp85,
              :combi_othGHG_N20_man_made_emissions_rcp3_u,
              Symbol("combi_othGHG_N20_man_made_emissions_rcp3_y[1]"),
              :othGHG_N20_man_made_emissions_rcp3,
              :combi_othGHG_N20_man_made_emissions_rcp45_u,
              Symbol("combi_othGHG_N20_man_made_emissions_rcp45_y[1]"),
              :othGHG_N20_man_made_emissions_rcp45,
              :combi_othGHG_N20_man_made_emissions_rcp6_u,
              Symbol("combi_othGHG_N20_man_made_emissions_rcp6_y[1]"),
              :othGHG_N20_man_made_emissions_rcp6,
              :combi_othGHG_N20_man_made_emissions_rcp85_u,
              Symbol("combi_othGHG_N20_man_made_emissions_rcp85_y[1]"),
              :othGHG_N20_man_made_emissions_rcp85,
              :combi_RCP_3_CO2_concentration_1850_2100_ppm_u,
              Symbol("combi_RCP_3_CO2_concentration_1850_2100_ppm_y[1]"),
              :RCP_3_CO2_concentration_1850_2100_ppm,
              :combi_RCP_45_CO2_concentration_1850_2100_ppm_u,
              Symbol("combi_RCP_45_CO2_concentration_1850_2100_ppm_y[1]"),
              :RCP_45_CO2_concentration_1850_2100_ppm,
              :combi_RCP_6_CO2_concentration_1850_2100_ppm_u,
              Symbol("combi_RCP_6_CO2_concentration_1850_2100_ppm_y[1]"),
              :RCP_6_CO2_concentration_1850_2100_ppm,
              :combi_RCP_85_CO2_concentration_1850_2100_ppm_u,
              Symbol("combi_RCP_85_CO2_concentration_1850_2100_ppm_y[1]"),
              :RCP_85_CO2_concentration_1850_2100_ppm,
              :combi_Montreal_gases_emissions_from_CO2e_C_Roads_u,
              Symbol("combi_Montreal_gases_emissions_from_CO2e_C_Roads_y[1]"),
            )
          end
          push!(variableConstructors, generateAlgebraicVariables4)
        end
        begin
          function generateAlgebraicVariables5()
            (
              :Montreal_gases_emissions_from_CO2e_C_Roads,
              :combi_Montreal_gases_emissions_from_CO2e_CAT_u,
              Symbol("combi_Montreal_gases_emissions_from_CO2e_CAT_y[1]"),
              :Montreal_gases_emissions_from_CO2e_CAT,
              :combi_Montreal_gases_emissions_pct_contribution_to_Total_CO2e_u,
              Symbol("combi_Montreal_gases_emissions_pct_contribution_to_Total_CO2e_y[1]"),
              :Montreal_gases_emissions_pct_contribution_to_Total_CO2e,
              :combi_N2O_man_made_emissions_from_CO2e_C_Roads_u,
              Symbol("combi_N2O_man_made_emissions_from_CO2e_C_Roads_y[1]"),
              :N2O_man_made_emissions_from_CO2e_C_Roads,
              :combi_N2O_man_made_emissions_from_CO2e_CAT_u,
              Symbol("combi_N2O_man_made_emissions_from_CO2e_CAT_y[1]"),
              :N2O_man_made_emissions_from_CO2e_CAT,
              :combi_N2O_emissions_pct_contribution_to_Total_CO2e_u,
              Symbol("combi_N2O_emissions_pct_contribution_to_Total_CO2e_y[1]"),
              :N2O_emissions_pct_contribution_to_Total_CO2e,
              :combi_Sea_level_rise_history_mm_u,
              Symbol("combi_Sea_level_rise_history_mm_y[1]"),
              :Sea_level_rise_history_mm,
              :combi_All_Human_activity_emissions_GtCO2e_yr_Base_for_tipping_point_search_u,
              Symbol("combi_All_Human_activity_emissions_GtCO2e_yr_Base_for_tipping_point_search_y[1]"),
              :combi_Arctic_freezing_cutoff_u,
              Symbol("combi_Arctic_freezing_cutoff_y[1]"),
              :combi_Blocked_by_H20_hist_Table_lookup_u,
              Symbol("combi_Blocked_by_H20_hist_Table_lookup_y[1]"),
              :combi_Blocked_by_H20_Table_lookup_u,
              Symbol("combi_Blocked_by_H20_Table_lookup_y[1]"),
              :combi_Effect_of_heat_in_atm_on_melting_ice__cut_off__u,
              Symbol("combi_Effect_of_heat_in_atm_on_melting_ice__cut_off__y[1]"),
              :combi_Exp_12a_reduction_in_emissions_LOOKUP_u,
              Symbol("combi_Exp_12a_reduction_in_emissions_LOOKUP_y[1]"),
              :combi_EXP_12b_CCS_from_2015_u,
              Symbol("combi_EXP_12b_CCS_from_2015_y[1]"),
              :combi_EXP_12e_white_surfaces_ease_in_u,
              Symbol("combi_EXP_12e_white_surfaces_ease_in_y[1]"),
              :combi_Fraction_blocked_by_CH4_spectrum_u,
              Symbol("combi_Fraction_blocked_by_CH4_spectrum_y[1]"),
              :combi_Fraction_blocked_by_CO2_spectrum_u,
              Symbol("combi_Fraction_blocked_by_CO2_spectrum_y[1]"),
              :combi_Future_shape_of_anthropogenic_aerosol_emissions_u,
              Symbol("combi_Future_shape_of_anthropogenic_aerosol_emissions_y[1]"),
              :combi_Melting_constraint_from_the_heat_in__ocean__surface_reservoir_u,
              Symbol("combi_Melting_constraint_from_the_heat_in__ocean__surface_reservoir_y[1]"),
              :combi_Melting_constraint_from_the_heat_in_atmosphere_reservoir_fraction_u,
              Symbol("combi_Melting_constraint_from_the_heat_in_atmosphere_reservoir_fraction_y[1]"),
              :combi_Melting_restraint_for_permafrost_from_heat_in_atmophere_u,
              Symbol("combi_Melting_restraint_for_permafrost_from_heat_in_atmophere_y[1]"),
              :combi_NATURE_CCS_removal_experiment_multiplier_u,
              Symbol("combi_NATURE_CCS_removal_experiment_multiplier_y[1]"),
              :combi_NF_clear_cut_fraction_u,
            )
          end
          push!(variableConstructors, generateAlgebraicVariables5)
        end
        begin
          function generateAlgebraicVariables6()
            (
              Symbol("combi_NF_clear_cut_fraction_y[1]"),
              :combi_NF_usage_cutoff_u,
              Symbol("combi_NF_usage_cutoff_y[1]"),
              :combi_Permafrost_melting_cutoff_u,
              Symbol("combi_Permafrost_melting_cutoff_y[1]"),
              :combi_RCPFossil_fuel_usage_cutoff_u,
              Symbol("combi_RCPFossil_fuel_usage_cutoff_y[1]"),
              :combi_Snowball_earth_cutoff_u,
              Symbol("combi_Snowball_earth_cutoff_y[1]"),
              :combi_Thermal_expansion_deep_in_1850_pct_u,
              Symbol("combi_Thermal_expansion_deep_in_1850_pct_y[1]"),
              :combi_Thermal_expansion_deep_pct_u,
              Symbol("combi_Thermal_expansion_deep_pct_y[1]"),
              :combi_Thermal_expansion_surface_in_1850_pct_u,
              Symbol("combi_Thermal_expansion_surface_in_1850_pct_y[1]"),
              :combi_Thermal_expansion_surface_pct_u,
              Symbol("combi_Thermal_expansion_surface_pct_y[1]"),
              :combi_TROP_deforestation_cutoff_u,
              Symbol("combi_TROP_deforestation_cutoff_y[1]"),
              :combi_TROP_deforestation_cutoff_effect_u,
              Symbol("combi_TROP_deforestation_cutoff_effect_y[1]"),
              :combi_TROP_deforestion_multiplier_wrt_2000_u,
              Symbol("combi_TROP_deforestion_multiplier_wrt_2000_y[1]"),
              :combi_Urbanzation_Effect_on_biomass_use_u,
              Symbol("combi_Urbanzation_Effect_on_biomass_use_y[1]"),
              :combi_Population_Lookup_bn_u,
              Symbol("combi_Population_Lookup_bn_y[1]"),
              :aux_1____Temp_gradient_minus_1___slope_,
              :Actual_time_to_degrade_all_GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_yr,
              :Actual_time_to_degrade_all_NF_Dead_biomass__litter_and_soil_organic_matter_SOM_yr,
              :Actual_time_to_degrade_all_TROP_Dead_biomass__litter_and_soil_organic_matter_SOM_yr,
              :Actual_time_to_degrade_all_TUNDRA_Dead_biomass__litter_and_soil_organic_matter_SOM_yr,
              :Aerosol_anthropogenic_emissions,
              :Albedo_Antartic,
              :Albedo_glacier,
              :Albedo_land_biomes,
              :Albedo_ocean_with_arctic_ice_changes,
              :Albedo_URBAN,
              :All_C_taken_out_due_to_change_in_land_use_GtC,
              :All_CH4_emissions_GtC_yr,
              :ALL_clouds_net_effect__pos_warming__neg_cooling__W_m2,
              :All_Human_activity_emissions_GtCO2e_yr_Base_for_tipping_point_search,
              :All_Human_activity_emissions_GtCO2e_yr,
              :All_N2O_emissions_MtN2O_yr,
              :Annual_flux_of_C_to_biomass_GtC_pr_yr,
              :Annual_glacial_ice_losing__pos__or_gaining__neg__GtIce_yr,
              :Annual_release_of_C_from_permafrost_GtC_y,
              :Antarctic_ice_area_decrease_Mkm2_pr_yr,
              :Antarctic_ice_area_increase_Mkm2_pr_yr,
              :Antarctic_ice_area_km2,
            )
          end
          push!(variableConstructors, generateAlgebraicVariables6)
        end
        begin
          function generateAlgebraicVariables7()
            (
              :Antarctic_ice_freezing_km3_yr,
              :Antarctic_ice_losing__pos__or_gaining__neg__GtIce_yr,
              :Antarctic_ice_melting__pos__or_freezing__neg__km3_yr,
              :Antarctic_ice_melting_as_water_km3_yr,
              :Antarctic_ice_melting_km3_yr,
              :Anthropogenic_aerosol_forcing,
              :Arctic_as_fraction_of_ocean,
              :Arctic_freezing_cutoff,
              :Arctic_ice_area_max_km2,
              :Arctic_ice_area_Mkm2,
              :Arctic_ice_melting__pos__or_freezing__neg__km2_yr,
              :Area_covered_by_high_clouds,
              :Area_covered_by_low_clouds,
              :Area_equivalent_of_linear_retreat_km2_yr,
              :Area_of_earth_Mkm2,
              :Area_of_land_Mkm2,
              :Atmos_heat_used_for_melting_1_yr,
              :Avg_C_concentration_in_top_layer,
              :Avg_CC_in_ocean_top_layer_ymoles_per_litre,
              :Avg_CO2_conc_in_ocean_top_layer_in_ppm,
              :Avg_earths_surface_albedo,
              :Avg_thickness_Antarctic_km,
              :Avg_thickness_glacier_km,
              :Avg_volcanic_activity_GtC_yr,
              :Barren_land_Mkm2,
              :BARREN_land_normal_albedo_Mkm2,
              :BARREN_land_white_Mkm2,
              :BB_radiation_at_atm_temp_in_atm_W_m2,
              :BB_radiation_at_surface_temp_ZJ_yr,
              :BB_radiation_at_Temp_in_atm_ZJ_yr,
              :Blocked_by_CH4,
              :Blocked_by_CO2,
              :Blocked_by_H20,
              :Blocked_by_H20_future_linear_equ,
              :Blocked_by_H20_future_poly_equ,
              :Blocked_by_H20_future_poly_equ_dyn,
              :Blocked_by_H20_future_poly_equ_dyn_0,
              :Blocked_by_H20_hist_Table_lookup,
              :Blocked_by_h20_poly_used,
              :Blocked_by_H20_Table_lookup,
              :Blocked_by_H2O_hist_and_fut,
              :Blocked_by_H2O_poly_dyn,
              :Blocked_by_H2O_poly_equ,
              :Blocked_by_otherGHG,
              :Blocking_multiplier_from_Kyoto_Flour,
              :Blocking_multiplier_from_Montreal_gases,
              :Blocking_multiplier_from_N2O,
              :Blocking_of_LW_rad_by_clouds,
              :C_absorption_by_ocean_from_atm_for_accumulation,
              :C_diffusion_into_ocean_from_atm,
            )
          end
          push!(variableConstructors, generateAlgebraicVariables7)
        end
        begin
          function generateAlgebraicVariables8()
            (
              :C_diffusion_into_ocean_from_atm_MtC_yr,
              :C_in_biomass,
              :C_in_GRASS_DeadB_and_soil_GtC,
              :C_in_GRASS_GtC,
              :C_in_GRASS_LB_GtC,
              :C_in_NF_DeadB_and_soil_GtC,
              :C_in_NF_GtC,
              :C_in_NF_LB_GtC,
              :C_in_TROP_DeadB_and_soil_GtC,
              :C_in_TROP_GtC,
              :C_in_TROP_LB_GtC,
              :C_in_TUNDRA_DeadB_and_soil_GtC,
              :C_in_TUNDRA_GtC,
              :C_in_TUNDRA_LB_GtC,
              :C_release_from_permafrost_melting_as_CO2_GtC_yr,
              :C_released_from_permafrost_as_either_CH4_or_CO2_in_GtC_yr,
              :C_removal_rate_from_atm_for_nature_May_2020_GtC_y,
              :C_runoff_from_biomass_soil,
              :Carbon_captured_and_stored_GtC___yr,
              :Carbon_concentration_in_cold_surface_ocean,
              :Carbon_concentration_in_CWTtB,
              :Carbon_concentration_in_deep_box_GtC_per_G_cubicM,
              :Carbon_concentration_in_intermdiate_box_GtC_per_G_cubicM,
              :Carbon_concentration_in_warm_surface,
              :Carbon_flow_from_cold_surface_downwelling_Gtc_per_yr,
              :Carbon_flow_from_cold_to_deep_GtC_per_yr,
              :Carbon_flow_from_deep,
              :Carbon_flow_from_intermediate_to_surface_box_GtC_per_yr,
              :Carbon_flow_from_warm_to_cold_surface_GtC_per_yr,
              :Carbon_in_cold_ocean_0_to_100m_1850_GtC,
              :Carbon_in_cold_ocean_trunk_100m_to_bottom_1850_GtC,
              :Carbon_in_ocean_deep_1k_to_bottom_ocean_1850_GtC,
              :Carbon_in_ocean_upwelling_100m_to_1km_1850_GtC,
              :Carbon_in_top_ocean_layer_1850_GtC,
              :Carbon_in_top_ocean_layer_GtC,
              :Carbon_in_warm_ocean_0_to_100m_1850_GtC,
              :CC_in_cold_downwelling_ymoles_per_litre,
              :CC_in_cold_downwelling_ymoles_per_litre__dimensionless_,
              :CC_in_cold_surface_ymoles_per_litre,
              :CC_in_cold_surface_ymoles_per_litre__dimensionless_,
              :CC_in_deep_box_ymoles_per_litre,
              :CC_in_deep_box_ymoles_per_litre__dimensionless_,
              :CC_in_intermediate_box_ymoles_per_litre,
              :CC_in_intermediate_box_ymoles_per_litre__dimensionless_,
              :CC_in_warm_surface_ymoles_per_litre,
              :CC_in_warm_surface_ymoles_per_litre__dimensionless_,
              :CH4_all_emissions_GtC_yr,
              :CH4_concentration_ppb,
              :CH4_conversion_to_CO2_and_H2O,
              :CH4_emissions_before_co2e_exp,
            )
          end
          push!(variableConstructors, generateAlgebraicVariables8)
        end
        begin
          function generateAlgebraicVariables9()
            (
              :CH4_emissions_CO2e_after_exp,
              :CH4_emissions_CO2e_after_exp_12a,
              :CH4_emissions_from_wetlands_destruction,
              :CH4_in_permafrost_area_melted___frozen_before_heat_constraint_GtC_yr,
              :CH4_in_the_atmosphere_converted_to_CO2,
              :CH4_per_sqkm_of_wetlands,
              :CH4_release___capture_from_permafrost_area_loss___gain_GtC_yr,
              :CO2_conc_atm_less_CO2_conc_sea,
              :CO2_conc_in_cold_surface_water_in_ppm,
              :CO2_conc_in_warm_surface_water_in_ppm,
              :CO2_concentration_calculated_as_a_1pct_pa_exponential_increase_ppm,
              :CO2_concentration_ppm,
              :CO2_concentration_used__after_any_experiments__ppm,
              :CO2_emissions_before_co2e_exp,
              :CO2_emissions_CO2e_after_exp,
              :CO2_flow_from_GRASS_to_atmosphere_GtC_yr,
              :CO2_flow_from_NF_to_atmosphere_GtC_yr,
              :CO2_flow_from_TROP_to_atmosphere_GtC_yr,
              :CO2_flow_from_TUNDRA_to_atmosphere_GtC_yr,
              :CO2_flux_from_atm_to_GRASS_for_new_growth_GtC_yr,
              :CO2_flux_from_atm_to_NF_for_new_growth_GtC_yr,
              :CO2_flux_from_atm_to_TROP_for_new_growth_GtC_yr,
              :CO2_flux_from_atm_to_TUNDRA_for_new_growth_GtC_yr,
              :CO2_flux_GRASS_to_atm_Gtc_yr,
              :CO2_flux_NF_to_atm_Gtc_yr,
              :CO2_flux_TROP_to_atm_GtC_yr,
              :CO2_flux_TUNDRA_to_atm_Gtc_yr,
              :CO2_radiative_forcing_since_1850_using_Myhre_formula_W_pr_m2,
              :Cold_dense_water_sinking_in_Sverdrup,
              :Concentration_of_C_in_ocean_top_layer_in_1850,
              :Contrib_of_BARREN_land_to_albedo_land,
              :Contrib_of_GRASS_to_albedo_land,
              :Contrib_of_ICE_ON_LAND_to_albedo_land,
              :Contrib_of_NF_to_albedo_land,
              :Contrib_of_TROP_to_albedo_land,
              :Contrib_of_TUNDRA_to_albedo_land,
              :Contribution_to_forcing_by_CH4,
              :Contribution_to_forcing_by_CO2,
              :Contribution_to_forcing_by_H2O,
              :Contribution_to_forcing_by_othGHG,
              :Convection_aka_sensible_heat_flow,
              :Convection_aka_sensible_heat_flow_W_m2,
              :Convection_as_f_of_temp_ZJ_yr,
              :Conversion_constant_GtC_to_ppm,
              :Conversion_constant_heat_ocean_deep_to_temp,
              :Conversion_heat_atm_to_temp,
              :Conversion_heat_surface_to_temp,
              :dbl_CO2_exp,
              :Deep_ocean__cold__volume,
              :delta_C_in_atmosphere_GtC,
            )
          end
          push!(variableConstructors, generateAlgebraicVariables9)
        end
        begin
          function generateAlgebraicVariables10()
            (
              :delta_C_in_biomass_GtC,
              :delta_C_in_ocean_GtC,
              :delta_CO2_concentration_since_1850_ppm,
              :delta_Temp_deep_ocean_degC,
              :Depositing_of_C_to_sediment,
              :Effect_of_acidification_on_NMPP,
              :Effect_of_C_concentration_on_NMPP,
              :Effect_of_CO2_on_new_biomass_growth,
              :Effect_of_heat_in_atm_on_melting_ice__cut_off_,
              :Effect_of_humidity_on_shifting_biomes,
              :Effect_of_population_and_urbanization_on_biomass_use,
              :Effect_of_temp_on_melting_antarctic_ice,
              :Effect_of_temp_on_melting_greenland_ice,
              :Effect_of_temp_on_melting_greenland_ice_that_slid_into_the_ocean,
              :Effect_of_temp_on_melting_or_freezing_glacial_ice,
              :Effect_of_temp_on_melting_or_freezing_of_Arctic_ice,
              :Effect_of_temperature_on_fire_incidence_dimensionless,
              :Effect_of_temperature_on_new_biomass_growth_dimensionless,
              :Effect_of_temperature_on_NMPP,
              :Effective_time_to_melt_Arctic_ice_at_the_reference_delta_temp,
              :Effective_time_to_melt_glacial_ice_at_the_reference_delta_temp,
              :Effective_time_to_melt_greenland_ice_at_the_reference_delta_temp,
              :Effective_time_to_melt_or_freeze_antarctic_ice_at_the_reference_delta_temp,
              :Effective_Time_to_regrow_TROP_after_deforesting_yr,
              :Emissions_of_aerosols_1850_to_2100_with_IPCC_Fig_pg_1037_Exp,
              :Emissions_of_anthro_CH4_1850_to_2100_GtC_yr,
              :Emissions_of_anthro_CH4_1850_to_2100_linearly_reduced_from_2015_GtC_yr,
              :Emissions_of_anthro_CH4_1850_to_2100_RCP_with_JR_shape_forecast_GtC_yr,
              :Emissions_of_anthro_CH4_1850_to_2100_RCP_with_Xpct_annual_incr_forecast_GtC_yr,
              :Emissions_of_anthro_CH4_1850_to_2100_RCP3_GtC_yr,
              :Emissions_of_anthro_CH4_1850_to_2100_RCP45_GtC_yr,
              :Emissions_of_anthro_CH4_1850_to_2100_RCP6_GtC_yr,
              :Emissions_of_anthro_CH4_1850_to_2100_RCP85_GtC_yr,
              :Emissions_of_anthro_CO2_1850_to_2100_linearly_reduced_from_2015_GtC_yr,
              :Emissions_of_anthro_CO2_1850_to_2100_RCP_with_Xpct_annual_incr_forecast_GtC_yr,
              :Emissions_of_CH4_1850_to_2100_GtC_yr_with_IPCC_Fig_pg_1037_Exp,
              :Emissions_of_CO2_1850_to_2100_GtC_yr_with_EXP_12a,
              :Emissions_of_CO2_1850_to_2100_GtC_yr_with_IPCC_Fig_pg_1037_Exp,
              :Emissions_of_CO2_1850_to_2100_GtC_yr,
              :Evaporation_aka_latent_heat_flow,
              :Evaporation_aka_latent_heat_flow_W_m2,
              :Evaporation_as_f_of_temp_ZJ_yr,
              :Exogenous_sliding_of_Greenland_ice_into_the_ocean,
              :Exp_12a_reduction_in_emissions,
              :Exp_12a_reduction_in_emissions_LOOKUP,
              :EXP_12b_CCS_from_2015,
              :EXP_12c_stopping_TROP_deforestation_from_2015,
              :EXP_12e_white_surfaces_ease_in,
              :exp0,
              :exp0_dyn,
            )
          end
          push!(variableConstructors, generateAlgebraicVariables10)
        end
        begin
          function generateAlgebraicVariables11()
            (
              :exp1,
              :exp1_dyn,
              :exp2,
              :exp2_dyn,
              :exp3,
              :exp3_dyn,
              :Experimental_doubling_of_constant_C_emissions,
              :Experimental_release_of_constant_fossil_C_emissions_GtC_yr,
              :Experimental_release_of_methane,
              :f_M_1750_N_2010__for_ch4_forcing,
              :f_M_2010_N_cur_,
              :f_M_cur_N_2010_,
              :f_M2010_N_1750__for_n20_forcing,
              :Flow_from_atm_to_biomass_GtC_pr_yr,
              :Flow_from_biomass_to_atm_Gtc_pr_yr,
              :Flow_of_cold_surface_water_welling_down_GcubicM_per_yr,
              :Flow_of_cold_water_sinking_to_very_bottom_GcubicM_per_yr,
              :Flow_of_heat_to_atm_ZJ_yr,
              :Flow_of_heat_to_deep_ocean,
              :Flow_of_heat_to_deep_ocean_btw_72_and_08,
              :Flow_of_heat_to_surface_ocean,
              :Flow_of_heat_to_surface_ocean_btw_1972_and_2008,
              :Flow_of_water_from_warm_to_cold_surface_Gcubicm_per_yr,
              :for_display_yr_on_yr_change_in_C_in_ocean_GtC_yr,
              :Frac_atm_absorption,
              :Frac_blocked_by_ALL_GHG,
              :Frac_blocked_by_ALL_GHG_LESS_watervapor,
              :Frac_vol_cold_ocean_0_to_100m_of_total,
              :Frac_vol_cold_ocean_downwelling_of_total,
              :Frac_vol_deep_ocean_of_total,
              :Frac_vol_ocean_upwelling_of_total,
              :Frac_vol_warm_ocean_0_to_100m_of_total,
              :Fraction_blocked_by_CH4_spectrum,
              :Fraction_blocked_by_CO2_spectrum,
              :Fraction_blocked_by_other_GHG,
              :Fraction_GRASS_being_deforested_1_yr,
              :Fraction_of_C_released_from_permafrost_released_as_CH4_dmnl,
              :Fraction_of_ocean_classified_as_cold_surface,
              :Fraction_TUNDRA_being_deforested_1_yr,
              :Future_shape_of_anthropogenic_aerosol_emissions,
              :Ga__BB_radiation_less_TOA_radiation_W_m2,
              :Glacial_ice_area_decrease_Mkm2_pr_yr,
              :Glacial_ice_area_increase_Mkm2_pr_yr,
              :Glacial_ice_area_km2,
              :Glacial_ice_freezing_km3_yr,
              :Glacial_ice_melting__pos__or_freezing__neg__km3_yr,
              :Glacial_ice_melting_as_water_km3_yr,
              :Glacial_ice_melting_km3_yr,
              :GRASS_being_deforested_Mkm2_yr,
              :GRASS_being_harvested_Mkm2_yr,
            )
          end
          push!(variableConstructors, generateAlgebraicVariables11)
        end
        begin
          function generateAlgebraicVariables12()
            (
              :GRASS_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr,
              :GRASS_Biomass_in_construction_material_being_burnt_GtBiomass_yr,
              :GRASS_Biomass_in_construction_material_left_to_rot_GtBiomass_yr,
              :GRASS_biomass_new_growing_GtBiomass___yr,
              :GRASS_burning_Mkm2_yr,
              :GRASS_Dead_biomass_decomposing_GtBiomass_yr,
              :GRASS_DeadB_and_SOM_tB_per_km2,
              :GRASS_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr,
              :GRASS_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr,
              :GRASS_for_construction_use_GtBiomass_yr,
              :GRASS_historical_deforestation_pct_yr,
              :GRASS_land_taken_out_of_use_GtBiomass,
              :GRASS_land_taken_out_of_use_Mkm2,
              :GRASS_living_biomass_densitiy_tBiomass_pr_km2,
              :GRASS_Living_biomass_rotting_GtBiomass_yr,
              :GRASS_potential_less_actual_living_biomass_GtBiomass,
              :GRASS_potential_living_biomass_GtBiomass,
              :GRASS_regrowing_after_being_burnt_Mkm2_yr,
              :GRASS_regrowing_after_being_deforested_Mkm2_yr,
              :GRASS_regrowing_after_harvesting_Mkm2_yr,
              :GRASS_runoff,
              :GRASS_soil_degradation_from_forest_fires_GtBiomass_yr,
              :GRASS_with_normal_cover_Mkm2,
              :Greenland_ice_area_decrease_Mkm2_pr_yr,
              :Greenland_ice_area_increase_Mkm2_pr_yr,
              :Greenland_ice_area_km2,
              :Greenland_ice_freezing_km3_yr,
              :Greenland_ice_losing__pos__or_gaining__neg__GtIce_yr,
              :Greenland_ice_melting__pos__or_freezing__neg__km3_yr,
              :Greenland_ice_melting_as_water_km3_yr,
              :Greenland_ice_melting_km3_yr,
              :Greenland_ice_melting_that_slid_into_the_ocean_km3_yr,
              :Greenland_ice_sliding_into_the_ocean_km3_yr,
              :Greenland_ice_that_slid_into_the_ocean_melting__pos__or_freezing__neg__km3_yr,
              :Guldberg_Waage_air_sea_formulation,
              :Heat_actually_gained___needed_for_freezing___unfreezing_of_permafrost_ZJ_yr,
              :Heat_flow_from_the_earths_core,
              :Heat_gained___needed_for_the_desired_freezing___unfreezing_of_permafrost_ZJ_yr,
              :Heat_used_in_melting__pos__or_freezing__neg__antarctic_ice_ZJ_yr,
              :Heat_used_in_melting__pos__or_freezing__neg__arctic_sea_ice_ZJ_yr,
              :Heat_used_in_melting__pos__or_freezing__neg__glacial_ice_ZJ_yr,
              :Heat_used_in_melting__pos__or_freezing__neg__Greenland_ice_that_slid_into_the_water_ZJ_yr,
              :Heat_used_in_melting__pos__or_freezing__neg__Greenland_ice_ZJ_yr,
              :Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_antarctic_ice_W_m2,
              :Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_antarctic_ice_ZJ_yr,
              :Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_arctic_ice_W_m2,
              :Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_arctic_ice_ZJ_yr,
              :Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_Greenland_ice_that_slid_into_the_ocean_W_m2,
              :Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_Greenland_ice_that_slid_into_the_ocean_ZJ_yr,
              :Heat_withdrawn_from_atm_by_melting__pos____added__neg__by_freezing_ice_W_m2,
            )
          end
          push!(variableConstructors, generateAlgebraicVariables12)
        end
        begin
          function generateAlgebraicVariables13()
            (
              :Heat_withdrawn_from_atm_by_melting__pos____added__neg__by_freezing_ice_ZJ_yr,
              :HI_clouds_net_effect__pos_warming__neg_cooling__W_m2,
              :Hist_Frac_atm_absorption,
              :Human_activity_CH4_emissions,
              :Human_activity_CH4_emissions_GtCO2e_yr,
              :Humidity_of_atmosphere_current_g_kg,
              :Humidity_of_atmosphere_g_kg,
              :Ice_on_land_area_Mkm2,
              :Incoming_solar_W_m2,
              :Incoming_solar_ZJ_yr,
              :InputEmissions_for_tipping_point_search,
              :Intercept_blocked_by_H20_future_equ,
              :Kyoto_Flour_concentration_ppt,
              :Kyoto_Flour_degradation,
              :Kyoto_Flour_emissions,
              :Kyoto_Flour_emissions_after_exp,
              :Kyoto_Flour_emissions_after_exp_12a,
              :Kyoto_Flour_emissions_before_exp,
              :Kyoto_Flour_emissions_GtCO2e_yr,
              :Kyoto_Flour_emissions_RCPs_or_JR52,
              :Land_area_km2,
              :Land_covered_with_ice_km2,
              :Land_covered_with_ice_Mkm2,
              :LO_clouds_net_effect__pos_warming__neg_cooling__W_m2,
              :LW_Blocking_multiplier_from_other_GHG,
              :LW_Clear_sky_emissions_from_atm,
              :LW_Clear_sky_emissions_from_atm_W_m2,
              :LW_clear_sky_emissions_to_surface,
              :LW_clear_sky_emissions_to_surface_W_m2,
              :LW_Cloudy_sky_emissions_from_atm,
              :LW_Cloudy_sky_emissions_from_atm_W_m2,
              :LW_HI_cloud_radiation,
              :LW_HI_cloud_radiation_reference_in_1850_W_m2,
              :LW_HI_cloud_radiation_W_m2,
              :LW_LO_cloud_radiation,
              :LW_LO_cloud_radiation_W_m2,
              :LW_radiation_blocked_by_CH4__pct_,
              :LW_radiation_blocked_by_CO2__pct_,
              :LW_radiation_blocked_by_H2O__pct_,
              :LW_radiation_blocked_by_other_GHG__pct_,
              :LW_re_radiated_by_clouds,
              :LW_re_radiated_by_clouds_W_m2,
              :LW_surface_emission,
              :LW_surface_emission_W_m2,
              :LW_surface_emissions_escaping_through_atm_window,
              :LW_surface_emissions_NOT_escaping_through_atm_window,
              :LW_surface_emissions_NOT_escaping_through_atm_window_W_m2,
              :LW_TOA_radiation_from_atm_to_space,
              :LW_TOA_radiation_from_atm_to_space_difference_wrt_1850,
              :LW_TOA_radiation_from_atm_to_space_difference_wrt_1850_W_m2,
            )
          end
          push!(variableConstructors, generateAlgebraicVariables13)
        end
        begin
          function generateAlgebraicVariables14()
            (
              :LW_TOA_radiation_from_atm_to_space_W_m2,
              :M_2010,
              :M_cur,
              :Man_made_CH4_emissions_pct,
              :Man_made_fossil_C_emissions_for_cumulation_GtC_yr,
              :Man_made_fossil_C_emissions_GtC_yr,
              :Man_made_fossil_C_emissions_GtCO2e_yr,
              :Melting_constraint_from_the_heat_in__ocean__surface_reservoir,
              :Melting_constraint_from_the_heat_in_atmosphere_reservoir_fraction,
              :Melting_restraint_for_permafrost_from_heat_in_atmophere,
              :Methane_hydrates_released_and_converted_to_CO2_by_bacteria_GtC,
              :Methanehydrate_experimental_release_GtC__yr,
              :MODEL_CH4_in_atm_in_ppb,
              :MODEL_CO2_concentration_in_atmosphere2_ppm,
              :Model_Volcanic_aerosol_forcing_W_m2,
              :Montreal_emissions_GtCO2e_yr,
              :Montreal_gases_concentration_ppt,
              :Montreal_gases_degradation,
              :Montreal_gases_emissions,
              :Montreal_gases_emissions_after_exp_12a,
              :Montreal_gases_emissions_before_exp,
              :Montreal_gases_emissions_CO2e_after_exp,
              :Montreal_gases_emissions_RCPs_or_JR52,
              :N_2010,
              :N_cur,
              :N20_emissions_RCPs_or_JR52,
              :N2O_concentration_ppb,
              :N2O_degradation_MtN2O_yr,
              :N2O_man_made_emissions,
              :N2O_man_made_emissions_after_exp,
              :N2O_man_made_emissions_exp_12a,
              :N2O_man_made_emissions_GtCO2e_yr,
              :NatEvent_d__slowing_down_ocean_circulation_from_2015,
              :Natural_CH4_emissions,
              :Natural_CH4_emissions_pct,
              :NATURE_CCS_Fig3_GtC_yr,
              :NATURE_CCS_removal_experiment_multiplier,
              :Net_additions_to_C_in_TUNDRA_DeadB_and_soil_GtC,
              :Net_additions_to_C_in_TUNDRA_LB_GtC,
              :Net_C_flow_from_atm_to_biomass_GtC_pr_yr,
              :Net_C_to_atm,
              :Net_C_to_atm_rate,
              :Net_CO2_flow_between_grass_and_atmosphere_GtC,
              :Net_CO2_flow_between_TUNDRA_and_atmosphere_GtC,
              :Net_flow_of_heat_into_surface,
              :Net_flux_to_ocean_GtC_yr,
              :Net_heat_flow_ocean_between_surface_and_deep_per_K_of_difference_ZJ_yr_K,
              :Net_heat_flow_ocean_from_surface_to_deep__ZJ_yr_,
              :Net_heat_flow_ocean_from_surface_to_deep_W_m2,
              :Net_heat_flow_to_atm_ZJ_yr__needed_for_comparisons_with_history_,
            )
          end
          push!(variableConstructors, generateAlgebraicVariables14)
        end
        begin
          function generateAlgebraicVariables15()
            (
              :Net_marine_primary_production_NMPP_GtC_pr_yr,
              :NEW_Temp_ocean_surface_in_1850_in_K,
              :NF_Avg_life_biomass_yr,
              :NF_being_deforested_Mkm2_yr,
              :NF_being_harvested_by_clear_cutting_Mkm2_yr,
              :NF_being_harvested_Mkm2_yr,
              :NF_being_harvested_normally_Mkm2_yr,
              :NF_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr,
              :NF_Biomass_in_construction_material_being_burnt_GtBiomass_yr,
              :NF_biomass_new_growing_GtBiomass___yr,
              :NF_burning_Mkm2_yr,
              :NF_clear_cut_fraction,
              :NF_Dead_biomass_decomposing_GtBiomass_yr,
              :NF_DeadB_and_SOM_tB_per_km2,
              :NF_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr,
              :NF_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr,
              :NF_for_construction_use_GtBiomass_yr,
              :NF_historical_deforestation_pct_yr,
              :NF_land_taken_out_of_use_GtBiomass,
              :NF_land_taken_out_of_use_Mkm2,
              :NF_living_biomass_densitiy_tBiomass_pr_km2,
              :NF_Living_biomass_rotting_GtBiomass_yr,
              :NF_potential_less_actual_living_biomass_GtBiomass,
              :NF_potential_living_biomass_GtBiomass,
              :NF_regrowing_after_being_burnt_Mkm2_yr,
              :NF_regrowing_after_being_clear_cut_Mkm2_yr,
              :NF_regrowing_after_being_deforested_Mkm2_yr,
              :NF_regrowing_after_harvesting_Mkm2_yr,
              :NF_runoff,
              :NF_soil_degradation_from_clear_cutting_GtBiomass_yr,
              :NF_soil_degradation_from_forest_fires_GtBiomass_yr,
              :NF_Speed_of_regrowth_yr,
              :NF_Sum_outflows_GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass_yr,
              :NF_TUNDRA_Biomass_in_construction_material_left_to_rot_GtBiomass_yr,
              :NF_usage_as_pct_of_potial_area,
              :NF_usage_cutoff,
              :NF_with_normal_cover_Mkm2,
              :Ocean_area_km2,
              :Ocean_circulation_slowdown_from_Greenland_ice_sliding_into_the_Atlantic,
              :Ocean_heat_used_for_melting_ZJ_yr,
              :Ocean_surface_area_km2,
              :Ocean_surface_delta_temp_to_1850_C,
              :Open_water_as_frac_of_ocean_area,
              :Outgoing_radiation_at_TOA_W_m2,
              :pct_change_in_fraction_blocked_by_ALL_GHG_wrt_1850,
              :pct_change_in_fraction_blocked_by_C02_wrt_1850,
              :pct_change_in_fraction_blocked_by_CH4_wrt_1850,
              :pct_change_in_fraction_blocked_by_othGHG_wrt_1850,
              :pct_reduction_in_C_in_GRASS,
              :pct_reduction_in_C_in_NF,
            )
          end
          push!(variableConstructors, generateAlgebraicVariables15)
        end
        begin
          function generateAlgebraicVariables16()
            (
              :pct_reduction_in_C_in_TROP,
              :pct_reduction_in_C_in_TUNDRA,
              :Permafrost_area_km2,
              :Permafrost_CH4_emissions_pct,
              :Permafrost_melting_cutoff,
              :pH_in_cold_deep_water,
              :ph_in_cold_downwelling_water,
              :pH_in_cold_suface_water,
              :pH_in_surface,
              :pH_in_upwelling_water,
              :pH_in_warm_surface_water,
              :POLICY_4_Stopping_logging_in_Northern_forests,
              :Radiation_balance_at_TOA_in_less_out_W_m2,
              :Radiative_forcing_from_CH4_wrt_1850_W_m2,
              :Radiative_forcing_from_CO2_wrt_1850_W_m2,
              :Radiative_forcing_from_H2O_wrt_1850_W_m2,
              :Radiative_forcing_from_othGHG_wrt_1850_W_m2,
              :Radiative_forcing_wrt_1850_W_m2_0,
              :Rate_of_destruction_of_wetlands,
              :Ratio_of_area_covered_by_high_clouds_current_to_1850,
              :Ratio_of_area_covered_by_low_clouds_current_to_1850,
              :RCPFossil_fuel_usage_cutoff,
              :Reflected_Solar_SW,
              :Reflected_Solar_SW_W_m2,
              :RF_CH4_IPCC_formula_W_m2,
              :RF_CO2_Model_Myhre_formula,
              :RF_CO2_Model_Myhre_formula_1850,
              :RF_CO2_RCP3_Myhre_formula,
              :RF_CO2_RCP45_Myhre_formula,
              :RF_CO2_RCP6_Myhre_formula,
              :RF_CO2_RCP85_Myhre_formula,
              :RF_from_Ga__BB_radiation_less_TOA_radiation_W_m2,
              :RF_N20_IPCC_formula_W_m2,
              :Sea_level_change_from_melting_ice_and_thermal_expansion_m,
              :Sea_level_change_from_thermal_expansion_deep_m,
              :Sea_level_change_from_thermal_expansion_surface_m,
              :Sea_level_rise_from_melting_ice_m,
              :Sea_level_rise_history_m,
              :Seconds_per_yr,
              :Sensitivity_of_high_cloud_coverage_to_temp,
              :Shifting_GRASS_to_DESERT_Mkm2_yr,
              :Shifting_GRASS_to_NF_Mkm2_yr,
              :Shifting_GRASS_to_TROP_Mkm2_yr,
              :Shifting_ice_on_land_to_tundra_Mkm2_yr,
              :Shifting_ice_to_tundra_from_detail_ice_on_land_Mkm2_pr_yr,
              :Shifting_NF_to_GRASS_Mkm2_yr,
              :Shifting_NF_to_TROP_Mkm2_yr,
              :Shifting_NF_to_Tundra_Mkm2_yr,
              :Shifting_TROP_to_GRASS_Mkm2_yr,
              :Shifting_TROP_to_NF_Mkm2_yr,
            )
          end
          push!(variableConstructors, generateAlgebraicVariables16)
        end
        begin
          function generateAlgebraicVariables17()
            (
              :Shifting_tundra_to_ice_from_detail_ice_on_land_Mkm2_pr_yr,
              :Shifting_tundra_to_ice_on_land_Mkm2_yr,
              :Shifting_Tundra_to_NF_Mkm2_yr,
              :SHUT_OFF_permafrost,
              :Sifting_DESERT_to_GRASS_Mkm2_yr,
              :Slider_for_H2O_slope,
              :Slope_blocked_by_H20_future_equ,
              :Slope_btw_temp_and_permafrost_melting___freezing,
              :Slope_of_effect_of_temp_shifting_DESERT_to_GRASS,
              :Slope_temp_vs_glacial_ice_melting,
              :Slowing_of_recapture_of_CH4_dmnl,
              :Snowball_earth_cutoff,
              :Solar_cycle_W_m2,
              :Solar_sine_forcing_W_m2,
              :Stop_of_human_deforestation,
              :Sum_biomes_Mkm2,
              :sum_blocked,
              :Sum_heat_to_ocean_1972_to_2008_ZJ,
              :Sum_outflows_GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass_yr,
              :Surface_deep__ocean__temp_diff_degC,
              :Surface_imbalance_pos_is_TO_surface,
              :Surface_imbalance_pos_is_TO_surface_W_m2,
              :Surface_ocean__warm__volume,
              :SW_Atmospheric_absorption,
              :SW_Atmospheric_absorption_W_m2,
              :SW_clear_sky_reflection_aka_scattering,
              :SW_clear_sky_reflection_aka_scattering_W_m2,
              :SW_HI_cloud_efffect_aka_cloud_albedo,
              :SW_HI_cloud_efffect_aka_TOA_albedo_W_m2,
              :SW_LO_cloud_efffect_aka_cloud_albedo,
              :SW_LO_cloud_efffect_aka_cloud_albedo_W_m2,
              :SW_surface_absorption,
              :SW_surface_absorption_W_m2_wrt_1850,
              :SW_surface_absorption_W_m2,
              :SW_surface_reflection,
              :SW_surface_reflection_W_m2_wrt_1850,
              :SW_surface_reflection_W_m2,
              :SW_to_surface,
              :SW_to_surface_W_m2,
              :Temp__ocean__deep_in_1850_in_K,
              :Temp__ocean__deep_in_C,
              :Temp__ocean__surface_in_K,
              :Temp_atm_average_K,
              :Temp_atm_in_C,
              :Temp_driver_to_shift_biomes_degC,
              :Temp_gradient,
              :Temp_gradient_minus_1,
              :Temp_gradient_minus_1___slope,
              :Temp_ocean_deep_in_K,
              :Temp_of_cold_downwelling_water,
            )
          end
          push!(variableConstructors, generateAlgebraicVariables17)
        end
        begin
          function generateAlgebraicVariables18()
            (
              :Temp_of_cold_surface_water,
              :Temp_surface_anomaly_compared_to_1850_degC,
              :Temp_surface_average_K,
              :Temp_surface_C,
              :Temp_surface_current_divided_by_value_in_1850_K_K,
              :Thermal_expansion_deep_in_1850_pct,
              :Thermal_expansion_deep_pct,
              :Thermal_expansion_surface_in_1850_pct,
              :Thermal_expansion_surface_pct,
              :Time_in_trunk,
              :Time_less_Greenland_slide_experiment_start_yr,
              :Time_to_degrade_Kyoto_Flour_yr,
              :Time_to_regrow_NF_after_buning_yr,
              :Tipping_point_search_emissions_GtCO2e_yr,
              :Tipping_point_year_of_peak,
              :Total_carbon_in_Ocean_1850_GtC,
              :Total_carbon_in_ocean_GtC,
              :Total_CO2e_emissions_as_f_peak__GtCO2e_yr,
              :Total_net_aerosol_forcing_ZJ_yr,
              :Total_net_aerosol_forcings_W_m2,
              :Total_sea_level_change_from_thermal_expansion_m,
              :Total_volume_of_ocean_water_GcubicM,
              :TROP_being_deforested_Mkm2_yr,
              :TROP_being_harvested_by_clear_cutting_Mkm2_yr,
              :TROP_being_harvested_Mkm2_yr,
              :TROP_being_harvested_normally_Mkm2_yr,
              :TROP_Biomass_in_construction_material_being_burnt_GtBiomass_yr,
              :TROP_Biomass_in_construction_material_left_to_rot_GtBiomass_yr,
              :TROP_biomass_new_growing_GtBiomass___yr,
              :TROP_burning_Mkm2_yr,
              :TROP_Dead_biomass_decomposing_GtBiomass_yr,
              :TROP_DeadB_and_SOM_tB_per_km2,
              :TROP_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr,
              :TROP_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr,
              :TROP_deforestation_cutoff,
              :TROP_deforestation_cutoff_effect,
              :TROP_deforested_as_pct_of_potial_area,
              :TROP_deforestion_multiplier_wrt_2000,
              :TROP_for_construction_use_GtBiomass_yr,
              :TROP_historical_deforestation_pct_yr,
              :TROP_land_taken_out_of_use_GtBiomass,
              :TROP_land_taken_out_of_use_Mkm2,
              :TROP_living_biomass_densitiy_tBiomass_pr_km2,
              :TROP_Living_biomass_rotting_GtBiomass_yr,
              :TROP_NF_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr,
              :TROP_NF_regrowing_after_being_burnt_Mkm2_yr,
              :TROP_NF_regrowing_after_harvesting_Mkm2_yr,
              :TROP_potential_less_actual_living_biomass_GtBiomass,
              :TROP_potential_living_biomass_GtBiomass,
              :TROP_regrowing_after_being_clear_cut_Mkm2_yr,
            )
          end
          push!(variableConstructors, generateAlgebraicVariables18)
        end
        begin
          function generateAlgebraicVariables19()
            (
              :TROP_regrowing_after_being_deforested_Mkm2_yr,
              :TROP_runoff,
              :TROP_runoff_time,
              :TROP_soil_degradation_from_clear_cutting_GtBiomass_yr,
              :TROP_soil_degradation_from_forest_fires_GtBiomass_yr,
              :TROP_Sum_outflows_GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass_yr,
              :TROP_Time_to_decompose_undisturbed_dead_biomass_yr,
              :TROP_Use_of_NF_biomass_for_energy_GtBiomass_yr,
              :TROP_with_normal_cover_Mkm2,
              :TUNDRA_being_deforested_Mkm2_yr,
              :TUNDRA_being_harvested_Mkm2_yr,
              :TUNDRA_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr,
              :TUNDRA_Biomass_in_construction_material_being_burnt_GtBiomass_yr,
              :TUNDRA_Biomass_in_construction_material_left_to_rot_GtBiomass_yr,
              :TUNDRA_biomass_new_growing_GtBiomass___yr,
              :TUNDRA_burning_Mkm2_yr,
              :TUNDRA_Dead_biomass_decomposing_GtBiomass_yr,
              :TUNDRA_DeadB_and_SOM_tB_per_km2,
              :TUNDRA_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr,
              :TUNDRA_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr,
              :TUNDRA_for_construction_use_GtBiomass_yr,
              :TUNDRA_historical_deforestation_pct_yr,
              :TUNDRA_land_taken_out_of_use_GtBiomass,
              :TUNDRA_land_taken_out_of_use_Mkm2,
              :TUNDRA_living_biomass_densitiy_tBiomass_pr_km2,
              :TUNDRA_Living_biomass_rotting_GtBiomass_yr,
              :TUNDRA_potential_less_actual_living_biomass_GtBiomass,
              :TUNDRA_potential_living_biomass_GtBiomass,
              :TUNDRA_regrowing_after_being_burnt_Mkm2_yr,
              :TUNDRA_regrowing_after_being_deforested_Mkm2_yr,
              :TUNDRA_regrowing_after_harvesting_Mkm2_yr,
              :TUNDRA_runoff,
              :TUNDRA_soil_degradation_from_forest_fires_GtBiomass_yr,
              :TUNDRA_Sum_outflows_GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass_yr,
              :TUNDRA_with_normal_cover_Mkm2,
              :UNIT_conversion_for_CH4_from_CO2e_to_C,
              :UNIT_conversion_for_CO2_from_CO2e_to_C,
              :UNIT_conversion_from_MtCH4_to_GtC,
              :UNIT_conversion_GtCO2e_to_GtC,
              :UNIT_conversion_mm_to_m,
              :UNIT_conversion_W_m2_earth_to_ZJ_yr,
              :UNIT_converter_GtC_Gm3_to_ymoles_litre,
              :Upper_to_deep_ocean_temp_diff_in_1850_degC,
              :Upwelling_from_deep,
              :Upwelling_to_surface,
              :Urban_area_fraction,
              :Urban_Mkm2,
              :Urbanzation_Effect_on_biomass_use,
              :Use_of_GRASS_biomass_for_construction_GtBiomass_yr,
              :Use_of_GRASS_biomass_for_energy_GtBiomass_yr,
            )
          end
          push!(variableConstructors, generateAlgebraicVariables19)
        end
        begin
          function generateAlgebraicVariables20()
            (
              :Use_of_GRASS_for_construction_in_2000_GtBiomass,
              :Use_of_GRASS_for_energy_in_2000_GtBiomass,
              :Use_of_NF_biomass_for_construction_GtBiomass_yr,
              :Use_of_NF_biomass_for_energy_GtBiomass_yr,
              :Use_of_NF_for_construction_in_2000_GtBiomass,
              :Use_of_NF_for_energy_in_2000_GtBiomass,
              :Use_of_TROP_biomass_for_construction_GtBiomass_yr,
              :Use_of_TROP_for_construction_in_2000_GtBiomass,
              :Use_of_TROP_for_energy_in_2000_GtBiomass,
              :Use_of_TUNDRA_biomass_for_construction_GtBiomass_yr,
              :Use_of_TUNDRA_biomass_for_energy_GtBiomass_yr,
              :Use_of_TUNDRA_for_construction_in_2000_GtBiomass,
              :Use_of_TUNDRA_for_energy_in_2000_GtBiomass,
              :Volcanic_aerosols_emissions,
              :Volcanic_aerosols_removed_from_stratosphere,
              :Volume_cold_ocean_0_to_100m,
              :Volume_cold_ocean_downwelling_100m_to_bottom,
              :Volume_expansion_from_thermal_expansion_deep_Gm3_km3,
              :Volume_expansion_from_thermal_expansion_surface_Gm3_km3,
              :Volume_ocean_deep_1km_to_bottom,
              :Volume_ocean_upwelling_100m_to_1km,
              :Volume_of_total_ocean_Gm3,
              :Volume_warm_ocean_0_to_100m,
              :Warming_due_to_CH4_blocking_W_m2,
              :Warming_due_to_CO2_blocking_W_m2,
              :Warming_due_to_othGHG_blocking_W_m2,
              :Warming_due_to_water_vapor_blocking_W_m2,
              :Years_of_exponential_rise_dless,
              :Years_of_exponential_rise_yr,
              :Years_still_needed_to_reach_zero_emission_goal_yr,
              :yr_on_yr_change_in_C_in_land_use_GtC_yr,
              :yr_on_yr_change_in_C_in_ocean_GtC_yr,
              :flow_NF_TUNDRA_Biomass_in_construction_material_left_to_rot_GtBiomass_yr,
              :flow_NF_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr,
              :flow_Evaporation_aka_latent_heat_flow,
              :flow_C_runoff_from_biomass_soil,
              :flow_Kyoto_Flour_degradation,
              :flow_N2O_degradation_MtN2O_yr,
              :flow_LW_TOA_radiation_from_atm_to_space,
              :flow_TROP_Living_biomass_rotting_GtBiomass_yr,
              :flow_CO2_flux_TUNDRA_to_atm_Gtc_yr,
              :flow_Sifting_DESERT_to_GRASS_Mkm2_yr,
              :flow_Upwelling_from_deep,
              :flow_TUNDRA_regrowing_after_being_burnt_Mkm2_yr,
              :flow_TUNDRA_runoff,
              :flow_Greenland_ice_that_slid_into_the_ocean_melting__pos__or_freezing__neg__km3_yr,
              :flow_NF_being_harvested_by_clear_cutting_Mkm2_yr,
              :flow_Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_arctic_ice_ZJ_yr,
              :flow_TROP_Biomass_in_construction_material_being_burnt_GtBiomass_yr,
              :flow_Glacial_ice_melting__pos__or_freezing__neg__km3_yr,
            )
          end
          push!(variableConstructors, generateAlgebraicVariables20)
        end
        begin
          function generateAlgebraicVariables21()
            (
              :flow_NATURE_CCS_Fig3_GtC_yr,
              :flow_NF_biomass_new_growing_GtBiomass___yr,
              :flow_LW_clear_sky_emissions_to_surface,
              :flow_CH4_in_the_atmosphere_converted_to_CO2,
              :flow_NF_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr,
              :flow_TROP_biomass_new_growing_GtBiomass___yr,
              :flow_GRASS_Living_biomass_rotting_GtBiomass_yr,
              :flow_TUNDRA_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr,
              :flow_TUNDRA_biomass_new_growing_GtBiomass___yr,
              :flow_CO2_flux_from_atm_to_NF_for_new_growth_GtC_yr,
              :flow_CH4_conversion_to_CO2_and_H2O,
              :flow_Flow_of_heat_to_deep_ocean_btw_72_and_08,
              :flow_GRASS_for_construction_use_GtBiomass_yr,
              :flow_Flow_of_cold_surface_water_welling_down_GcubicM_per_yr,
              :flow_TROP_for_construction_use_GtBiomass_yr,
              :flow_Annual_glacial_ice_losing__pos__or_gaining__neg__GtIce_yr,
              :flow_NF_runoff,
              :flow_NF_soil_degradation_from_forest_fires_GtBiomass_yr,
              :flow_GRASS_runoff,
              :flow_Greenland_ice_sliding_into_the_ocean_km3_yr,
              :flow_TUNDRA_soil_degradation_from_forest_fires_GtBiomass_yr,
              :flow_Greenland_ice_melting__pos__or_freezing__neg__km3_yr,
              :flow_SW_surface_absorption,
              :flow_All_N2O_emissions_MtN2O_yr,
              :flow_NF_being_harvested_normally_Mkm2_yr,
              :flow_Kyoto_Flour_emissions,
              :flow_CO2_flux_from_atm_to_TUNDRA_for_new_growth_GtC_yr,
              :flow_Shifting_NF_to_TROP_Mkm2_yr,
              :flow_GRASS_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr,
              :flow_TUNDRA_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr,
              :flow_Carbon_flow_from_warm_to_cold_surface_GtC_per_yr,
              :flow_Shifting_GRASS_to_DESERT_Mkm2_yr,
              :flow_NF_being_deforested_Mkm2_yr,
              :flow_Heat_withdrawn_from_atm_by_melting__pos____added__neg__by_freezing_ice_ZJ_yr,
              :flow_GRASS_biomass_new_growing_GtBiomass___yr,
              :flow_Man_made_fossil_C_emissions_GtC_yr,
              :flow_Greenland_ice_melting_as_water_km3_yr,
              :flow_TROP_runoff,
              :flow_Antarctic_ice_losing__pos__or_gaining__neg__GtIce_yr,
              :flow_Carbon_flow_from_cold_surface_downwelling_Gtc_per_yr,
              :flow_NF_regrowing_after_harvesting_Mkm2_yr,
              :flow_TROP_Dead_biomass_decomposing_GtBiomass_yr,
              :flow_TUNDRA_being_deforested_Mkm2_yr,
              :flow_Shifting_TROP_to_GRASS_Mkm2_yr,
              :flow_CH4_release___capture_from_permafrost_area_loss___gain_GtC_yr,
              :flow_Volcanic_aerosols_emissions,
              :flow_GRASS_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr,
              :flow_Natural_CH4_emissions,
              :flow_Flow_of_heat_to_atm_ZJ_yr,
              :flow_NF_Biomass_in_construction_material_being_burnt_GtBiomass_yr,
            )
          end
          push!(variableConstructors, generateAlgebraicVariables21)
        end
        begin
          function generateAlgebraicVariables22()
            (
              :flow_Flow_of_heat_to_deep_ocean,
              :flow_LW_surface_emission,
              :flow_NF_regrowing_after_being_burnt_Mkm2_yr,
              :flow_TROP_NF_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr,
              :flow_TUNDRA_Biomass_in_construction_material_being_burnt_GtBiomass_yr,
              :flow_Man_made_fossil_C_emissions_for_cumulation_GtC_yr,
              :flow_C_absorption_by_ocean_from_atm_for_accumulation,
              :flow_Antarctic_ice_melting__pos__or_freezing__neg__km3_yr,
              :flow_Annual_flux_of_C_to_biomass_GtC_pr_yr,
              :flow_GRASS_Biomass_in_construction_material_being_burnt_GtBiomass_yr,
              :flow_NF_regrowing_after_being_deforested_Mkm2_yr,
              :flow_Greenland_ice_losing__pos__or_gaining__neg__GtIce_yr,
              :flow_CO2_flux_from_atm_to_TROP_for_new_growth_GtC_yr,
              :flow_NF_soil_degradation_from_clear_cutting_GtBiomass_yr,
              :flow_Annual_release_of_C_from_permafrost_GtC_y,
              :flow_Avg_volcanic_activity_GtC_yr,
              :flow_TUNDRA_regrowing_after_harvesting_Mkm2_yr,
              :flow_Shifting_ice_on_land_to_tundra_Mkm2_yr,
              :flow_C_diffusion_into_ocean_from_atm,
              :flow_Glacial_ice_melting_as_water_km3_yr,
              :flow_NF_for_construction_use_GtBiomass_yr,
              :flow_Flow_of_heat_to_surface_ocean,
              :flow_TROP_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr,
              :flow_C_released_from_permafrost_as_either_CH4_or_CO2_in_GtC_yr,
              :flow_TROP_soil_degradation_from_forest_fires_GtBiomass_yr,
              :flow_TROP_being_harvested_by_clear_cutting_Mkm2_yr,
              :flow_NF_regrowing_after_being_clear_cut_Mkm2_yr,
              :flow_GRASS_being_harvested_Mkm2_yr,
              :flow_Convection_aka_sensible_heat_flow,
              :flow_TUNDRA_for_construction_use_GtBiomass_yr,
              :flow_NF_burning_Mkm2_yr,
              :flow_NF_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr,
              :flow_TUNDRA_burning_Mkm2_yr,
              :flow_CO2_flux_TROP_to_atm_GtC_yr,
              :flow_Shifting_tundra_to_ice_on_land_Mkm2_yr,
              :flow_Flow_of_water_from_warm_to_cold_surface_Gcubicm_per_yr,
              :flow_Shifting_Tundra_to_NF_Mkm2_yr,
              :flow_Flow_of_heat_to_surface_ocean_btw_1972_and_2008,
              :flow_TUNDRA_Living_biomass_rotting_GtBiomass_yr,
              :flow_Methanehydrate_experimental_release_GtC__yr,
              :flow_GRASS_regrowing_after_being_burnt_Mkm2_yr,
              :flow_Montreal_gases_degradation,
              :flow_Carbon_flow_from_cold_to_deep_GtC_per_yr,
              :flow_GRASS_soil_degradation_from_forest_fires_GtBiomass_yr,
              :flow_Shifting_TROP_to_NF_Mkm2_yr,
              :flow_GRASS_being_deforested_Mkm2_yr,
              :flow_Shifting_GRASS_to_NF_Mkm2_yr,
              :flow_TROP_being_deforested_Mkm2_yr,
              :flow_Arctic_ice_melting__pos__or_freezing__neg__km2_yr,
              :flow_CO2_flux_from_atm_to_GRASS_for_new_growth_GtC_yr,
            )
          end
          push!(variableConstructors, generateAlgebraicVariables22)
        end
        begin
          function generateAlgebraicVariables23()
            (
              :flow_GRASS_regrowing_after_being_deforested_Mkm2_yr,
              :flow_Net_C_to_atm_rate,
              :flow_Methane_hydrates_released_and_converted_to_CO2_by_bacteria_GtC,
              :flow_LW_surface_emissions_NOT_escaping_through_atm_window,
              :flow_Antarctic_ice_melting_as_water_km3_yr,
              :flow_TUNDRA_Biomass_in_construction_material_left_to_rot_GtBiomass_yr,
              :flow_TROP_NF_regrowing_after_harvesting_Mkm2_yr,
              :flow_TUNDRA_being_harvested_Mkm2_yr,
              :flow_Net_heat_flow_ocean_from_surface_to_deep__ZJ_yr_,
              :flow_TROP_regrowing_after_being_clear_cut_Mkm2_yr,
              :flow_GRASS_Biomass_in_construction_material_left_to_rot_GtBiomass_yr,
              :flow_TROP_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr,
              :flow_Carbon_flow_from_deep,
              :flow_Rate_of_destruction_of_wetlands,
              :flow_Montreal_gases_emissions,
              :flow_LW_re_radiated_by_clouds,
              :flow_Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_antarctic_ice_ZJ_yr,
              :flow_Depositing_of_C_to_sediment,
              :flow_TUNDRA_Dead_biomass_decomposing_GtBiomass_yr,
              :flow_TUNDRA_regrowing_after_being_deforested_Mkm2_yr,
              :flow_TROP_burning_Mkm2_yr,
              :flow_TROP_NF_regrowing_after_being_burnt_Mkm2_yr,
              :flow_SW_Atmospheric_absorption,
              :flow_GRASS_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr,
              :flow_GRASS_regrowing_after_harvesting_Mkm2_yr,
              :flow_TROP_being_harvested_normally_Mkm2_yr,
              :flow_C_release_from_permafrost_melting_as_CO2_GtC_yr,
              :flow_Human_activity_CH4_emissions,
              :flow_GRASS_Dead_biomass_decomposing_GtBiomass_yr,
              :flow_TROP_Biomass_in_construction_material_left_to_rot_GtBiomass_yr,
              :flow_TROP_soil_degradation_from_clear_cutting_GtBiomass_yr,
              :flow_TUNDRA_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr,
              :flow_Shifting_NF_to_GRASS_Mkm2_yr,
              :flow_Heat_flow_from_the_earths_core,
              :flow_Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_Greenland_ice_that_slid_into_the_ocean_ZJ_yr,
              :flow_TROP_regrowing_after_being_deforested_Mkm2_yr,
              :flow_C_removal_rate_from_atm_for_nature_May_2020_GtC_y,
              :flow_GRASS_burning_Mkm2_yr,
              :flow_CO2_flux_GRASS_to_atm_Gtc_yr,
              :flow_Upwelling_to_surface,
              :flow_NF_Dead_biomass_decomposing_GtBiomass_yr,
              :flow_Carbon_captured_and_stored_GtC___yr,
              :flow_Volcanic_aerosols_removed_from_stratosphere,
              :flow_Carbon_flow_from_intermediate_to_surface_box_GtC_per_yr,
              :flow_Greenland_ice_melting_that_slid_into_the_ocean_km3_yr,
              :flow_Shifting_NF_to_Tundra_Mkm2_yr,
              :flow_Shifting_GRASS_to_TROP_Mkm2_yr,
              :flow_NF_Living_biomass_rotting_GtBiomass_yr,
              :flow_CO2_flux_NF_to_atm_Gtc_yr,
              :flow_Flow_of_cold_water_sinking_to_very_bottom_GcubicM_per_yr,
            )
          end
          push!(variableConstructors, generateAlgebraicVariables23)
        end
        begin
          function generateAlgebraicVariables24()
            (
              :flow_Biological_removal_of_C_from_WSW_GtC_per_yr,
              :C_in_ocean_1_yr_ago_GtC_DL,
              :Atmos_heat_used_for_melting_last_year_1_yr,
              :Ocean_heat_used_for_melting_last_year_ZJ_yr,
              :C_in_atm_1_yr_ago_GtC,
              :C_in_atm_1_yr_ago_GtC_RT1,
              :C_in_atm_1_yr_ago_GtC_RT2,
              :C_in_atm_1_yr_ago_GtC_DL,
              :All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC,
              :All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_RT1,
              :All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_RT2,
              :All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_DL,
              :ifEq_tmp304,
              :ifEq_tmp305,
            )
          end
          push!(variableConstructors, generateAlgebraicVariables24)
        end
      end
      allVariables = []
      for constructor in variableConstructors
        t = Symbolics.variable(:t, T = Real)
        vars = map((n -> (n, (Symbolics.variable(n, T = Symbolics.FnType{Tuple{Real},Real}))(t))), constructor())
        push!(allVariables, vars)
      end
      vars = collect(Iterators.flatten(allVariables))
      for (sym, var) in vars
        eval(:($sym = $var))
      end
      local irreductableSyms = [:Time, :Kyoto_Flour_concentration_ppt, :ifEq_tmp304, :ifEq_tmp304, :Time, :Montreal_gases_concentration_ppt, :ifEq_tmp305, :ifEq_tmp305]
      for sym in irreductableSyms
        eval(:($sym = SymbolicUtils.setmetadata($sym, ModelingToolkit.VariableIrreducible, true)))
      end
      vars = map((x -> last(x)), vars)
      pars = Dict(
        Future_volcanic_emissions => 0.0,
        Albedo_Antarctic_hist => 0.7,
        Albedo_Antarctic_sens => 0.7,
        Albedo_BARREN_normal => 0.17,
        Albedo_BARREN_white => 0.7,
        Albedo_DESERT_normal => 0.24,
        Albedo_glacier_hist => 0.4,
        Albedo_glacier_sens => 0.4,
        Albedo_GRASS_burnt => 0.08,
        Albedo_GRASS_deforested => 0.3,
        Albedo_GRASS_normal_cover => 0.16,
        Albedo_Greenland => 0.7,
        Albedo_NF_burnt => 0.13,
        Albedo_NF_deforested => 0.18,
        Albedo_NF_normal_cover => 0.08,
        Albedo_TROP_burnt => 0.1,
        Albedo_TROP_deforested => 0.168,
        Albedo_TROP_normal_cover => 0.14,
        Albedo_TUNDRA_burnt => 0.23,
        Albedo_TUNDRA_deforested => 0.23,
        Albedo_TUNDRA_normal_cover => 0.23,
        Albedo_URBAN_normal => 0.15,
        Amount_methane_hydrates__clathrates__experimentally_released_GtC => 0.0,
        Amt_of_constant_emissions_GtC_yr => 4.0,
        Annual_pct_increase_CH4_emissions_from_2015_pct_yr => 0.0,
        Annual_pct_increase_CO2_emissions_from_2015_pct_yr => 0.0,
        Antarctic_ice_volume_in_1850_km3 => 3.0e7,
        Arctic_ice_albedo_1850 => 0.7,
        Arctic_ice_area_in_1850_km2 => 1.34e7,
        Arctic_surface_temp_delay_yr => 15.0,
        Area_covered_by_high_clouds_in_1850 => 0.2,
        Area_covered_by_low_clouds_in_1850 => 0.4,
        Area_equivalent_of_1km_linear_retreat_km2 => 17500.0,
        Area_of_earth_m2 => 5.1e14,
        Area_of_ocean_at_surface_361900_Gm2 => 361900.0,
        Atmos_heat_used_for_melting_Initially_1_yr => 0.0,
        Average_thickness_arctic_ice_km => 0.0025,
        Avg_amount_of_C_in_the_form_of_CH4_per_km2 => 4.8e-5,
        Avg_depth_of_permafrost_km => 0.1,
        Avg_flatness_of_worlds_coastline => 1.0,
        Avg_thickness_Antarctic_hist_km => 2.14,
        Avg_thickness_Antarctic_sens_km => 2.14,
        Avg_thickness_Greenland_km => 1.35,
        C_in_atmosphere_in_1850_GtC => 600.0,
        C_in_the_form_of_CH4_in_atm_1850 => 1.69,
        Carbon_per_biomass_tC_per_tBiomass => 0.5,
        CC_in_cold_ocean_0_to_100m_1850_ymoles_per_litre => 2240.0,
        CC_in_cold_ocean_downwelling_100m_bottom_1850_ymoles_per_litre => 2240.0,
        CC_in_ocean_upwelling_100m_to_1km_1850_ymoles_per_litre => 2240.0,
        CC_in_warm_ocean_0_to_100m_1850_ymoles_per_litre => 2240.0,
        CC_ocean_deep_1km_to_bottom_1850_ymoles_per_litre => 2240.0,
        CH4_concentration_in_2010_ppb => 1720.81,
        CH4_halflife_in_atmosphere => 7.3,
        Cold_dense_water_sinking_in_Sverdrup_in_1850 => 35.0,
        Constant_anthropogenic_CH4_emissions => 0.2,
        Convection_as_f_of_incoming_solar_in_1850 => 0.071,
        conversion_factor_CH4_Gt_to_ppb => 468.0,
        Conversion_from_Kyoto_Flour_amount_to_concentration_ppt_kt => 0.04,
        Conversion_from_Montreal_gases_amount_to_concentration_ppt_kt => 0.04,
        Conversion_Millionkm2_to_km2_Mkm2_km2 => 1.0e-6,
        Conversion_of_anthro_aerosol_emissions_to_forcing => -1.325,
        Conversion_of_volcanic_aerosol_emissions_to_CO2_emissions_GtC_pr_VAE => 2.8,
        Conversion_of_volcanic_aerosol_forcing_to_volcanic_aerosol_emissions => -1.0,
        Conversion_ymoles_per_kg_to_pCO2_yatm => 0.127044,
        Densitiy_of_water_relative_to_ice => 0.916,
        Duration_of_destruction_yr => 5.0,
        Emissions_of_natural_CH4_GtC_yr => 0.19,
        Emissivity_atm => 1.0,
        Emissivity_surface => 1.0,
        Evaporation_as_fraction_of_incoming_solar_in_1850 => 0.289,
        EXP_12f_Stratospheric_scattering_experiment_0_off_1_on => float(0),
        Experimental_doubling_of_constant_C_emissions_how_long_yr => 5.0,
        Experimental_doubling_of_constant_C_emissions_how_much_1_100pct => 0.0,
        Experimental_doubling_of_constant_C_emissions_when_yr => 30000.0,
        Frac_of_surface_emission_through_atm_window => 0.051,
        Frac_SW_clear_sky_reflection_aka_scattering => 0.0837,
        Frac_SW_HI_cloud_efffect_aka_cloud_albedo => 0.006,
        Frac_SW_LO_cloud_efffect_aka_cloud_albedo => 0.158,
        Fraction_of_C_released_from_permafrost_released_as_CH4_hist_dmnl => 1.0,
        Fraction_of_C_released_from_permafrost_released_as_CH4_sensitivity_dmnl => 1.0,
        Fraction_of_earth_surface_as_ocean => 0.7,
        Fraction_of_heat_needed_to_melt_antarctic_ice_coming_from_air => 0.6,
        Fraction_of_heat_needed_to_melt_arctic_ice_coming_from_air => 0.5,
        Fraction_of_heat_needed_to_melt_Greenland_ice_that_slid_into_the_ocean_coming_from_air => 0.1,
        Fraction_of_methane_hydrates_released_from_the_ocean_converted_to_CO2_before_it_is_relased_to_the_atmosphere => 0.9,
        Fraction_of_ocean_classified_warm_surface => 0.8,
        Glacial_ice_volume_in_1850_km3 => 167000.0,
        Global_Warming_Potential_CH4 => 25.0,
        Global_Warming_Potential_N20 => 298.0,
        GRASS_area_burned_in_1850_Mkm2 => 1.0,
        GRASS_area_deforested_in_1850_Mkm2 => 0.5,
        GRASS_area_harvested_in_1850_Mkm2 => 2.5,
        GRASS_Avg_life_biomass_yr => 100.0,
        GRASS_Avg_life_of_building_yr => 10.0,
        GRASS_Biomass_locked_in_construction_material_in_1850_GtBiomass => 1.5,
        GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_in_1850_GtBiomass => 1200.0,
        GRASS_Fraction_of_construction_waste_burned_0_1 => 0.5,
        GRASS_fraction_of_DeadB_and_SOM_being_destroyed_by_deforesting => 1.0,
        GRASS_fraction_of_DeadB_and_SOM_being_destroyed_by_energy_harvesting => 0.1,
        GRASS_fraction_of_DeadB_and_SOM_destroyed_by_natural_fires => 0.0,
        GRASS_living_biomass_densitiy_in_1850_tBiomass_pr_km2 => 14500.0,
        GRASS_Living_biomass_in_1850_GtBiomass => 310.0,
        GRASS_Normal_fire_incidence_fraction_yr => 1.0,
        GRASS_Ref_historical_deforestation_pct_yr => 0.1,
        GRASS_runoff_time => 2000.0,
        GRASS_Speed_of_regrowth_yr => 2.0,
        GRASS_Time_to_decompose_undisturbed_dead_biomass_yr => 1000.0,
        Greenland_ice_slide_circulation_slowdown_effect => 0.33,
        Greenland_ice_volume_in_1850_km3 => 2.93e6,
        Greenland_slide_experiment_how_much_sildes_in_the_ocean_fraction => 0.25,
        Greenland_slide_experiment_over_how_many_years_yr => 70.0,
        GtIce_vs_km3 => 0.9167,
        Heat_gained___needed_to_freeze___unfreeze_1_km3_permafrost_ZJ_km3 => 0.0001717,
        Heat_in__ocean__deep_in_1850_ZJ => 1.9532e6,
        Heat_in_atmosphere_in_1850_ZJ => 1025.67,
        Heat_in_surface_in_1850_ZJ => 25000.0,
        Heat_needed_to_melt_1_km3_of_ice_ZJ => 0.0003327,
        Hist_Avg_thickness_glacier_km => 0.23,
        Hist_Net_heat_flow_ocean_between_surface_and_deep_per_K_of_difference_ZJ_yr_K => 10.0,
        Hist_NF_Avg_life_biomass_yr => 60.0,
        Hist_NF_Speed_of_regrowth_yr => 3.0,
        Hist_Slope_of_effect_of_temp_shifting_DESERT_to_GRASS => 0.4,
        Hist_Slope_temp_vs_glacial_ice_melting => 1.0,
        Hist_Time_in_trunk => 234.638,
        Hist_Time_to_degrade_Kyoto_Flour_yr => 50.0,
        Hist_Time_to_regrow_NF_after_buning_yr => 30.0,
        Hist_TROP_runoff_time => 2000.0,
        Hist_TROP_Time_to_decompose_undisturbed_dead_biomass_yr => 24.0,
        K_to_C_conversion_C_K => 273.15,
        Kyoto_Flour_Global_Warming_Potential => 7000.0,
        Land_surface_temp_adjustment_time_yr => 25.0,
        LW_ALL_cloud_radiation_reference_in_1850_W_m2 => 27.9,
        LW_LO_cloud_radiation_reference_in_1850_W_m2 => 20.0,
        LW_radiation_fraction_blocked_by_other_GHG_in_1850 => 0.0398,
        Man_made_CH4_emissions_in_2015_GtC => 0.303,
        Man_made_CO2_emissions_in_2015_GtC => 10.0,
        MAX_NATURE_CCS_removal_in_2050_GtCO2e_yr => 35.0,
        Melting_of_permafrost_at_all_depths_at_4_deg_C_temp_diff_km_yr => 0.71,
        Montreal_Global_Warming_Potential => 10000.0,
        Myhre_constant_for_CH4 => 0.0594,
        Myhre_constant_for_CO2 => 5.35,
        Myhre_constant_for_N20 => 0.12,
        N2O_concentration_in_2010_ppb => 363.504,
        N2O_in_atmosphere_MtN2O_in_1850 => 900.0,
        N2O_natural_emissions => 9.0,
        Net_marine_primary_production_in_1850 => 0.4,
        NEvt_13a_double_rate_of_melting_ice_and_permafrost => float(1),
        NEvt_13b2_Double_incidence_of_biomass_fires => float(1),
        NEvt_13b3_double_sunspot_amplitude_from_2015_onwards_1_normal_2_double => float(1),
        NEvt_13c1_increase_in_area_covered_by_low_clouds => float(1),
        NEvt_13d_Greenland_slide_experiment_start_yr => float(3000000),
        NEvt_2a_Volcanic_eruptions_in_the_future_VAEs_first_future_pulse => float(21000000),
        NEvt_3b_increase_in_area_covered_by_high_clouds => float(1),
        NF_area_burned_in_1850_Mkm2 => 2.5,
        NF_area_deforested_in_1850_Mkm2 => 0.0,
        NF_area_harvested_in_1850_Mkm2 => 1.0,
        NF_Avg_life_of_building_yr => 20.0,
        NF_Biomass_locked_in_construction_material_in_1850_GtBiomass => 3.0,
        NF_Dead_biomass__litter_and_soil_organic_matter_SOM_in_1850_GtBiomass => 330.0,
        NF_DeadB_and_SOM_densitiy_in_1850_tBiomass_pr_km2 => 27500.0,
        NF_Fraction_of_construction_waste_burned_0_1 => 0.5,
        NF_fraction_of_DeadB_and_SOM_being_destroyed_by_clear_cutting => 0.5,
        NF_fraction_of_DeadB_and_SOM_being_destroyed_by_deforesting => 1.0,
        NF_fraction_of_DeadB_and_SOM_being_destroyed_by_energy_harvesting => 0.1,
        NF_fraction_of_DeadB_and_SOM_destroyed_by_natural_fires => 0.0,
        NF_living_biomass_densitiy_in_1850_tBiomass_pr_km2 => 7500.0,
        NF_Living_biomass_in_1850_GtBiomass => 115.0,
        NF_Normal_fire_incidence_fraction_yr => 0.7,
        NF_Ref_historical_deforestation_pct_yr => 0.02,
        NF_runoff_time => 2000.0,
        NF_Time_to_decompose_undisturbed_dead_biomass_yr => 250.0,
        Ocean_heat_used_for_melting_Initially_1_yr => 0.0,
        Ocean_slowdown_experimental_factor => 1.0,
        Open_ocean_albedo => 0.065,
        Over_how_many_yrs_methane_hydrate_release_yr => 5.0,
        per_annum_yr => 1.0,
        Policy_1_Reducing_GHG_emissions_by_one_third_by_2035 => float(0),
        Policy_2_Large_scale_implementation_of_carbon_capture_and_geological_storage__CCS_ => float(0),
        Population_2000_bn => 6.1,
        Pressure_adjustment_deep_pct => 1.0,
        Pressure_adjustment_surface_pct => 0.2,
        Rate_of_wetland_destruction_pct_of_existing_wetlands_yr => 0.0,
        Ratio_of_methane_in_tundra_to_wetland => 4.0,
        Ref_shifting_biome_yr => 50.0,
        Ref_temp_difference__4_degC_ => 4.0,
        Ref_temp_difference_for_antarctic_ice_melting__3_degC_ => 3.0,
        Ref_temp_difference_for_Arctic_ice_melting => 0.4,
        Ref_temp_difference_for_glacial_ice_melting__1_degC_ => 3.0,
        Ref_temp_difference_for_greenland_ice_melting_C => 1.0,
        Ref_temp_difference_for_greenland_ice_that_slid_into_the_ocean_melting_degC => 1.0,
        Reference_temp_C => 10.0,
        Reference_Time_to_regrow_TROP_after_deforesting_yr => 10000.0,
        SCALE_and_UNIT_converter_zero_C_to_K => 273.15,
        Sens_Avg_thickness_glacier_km => 0.23,
        Sens_Frac_atm_absorption => 0.220588,
        Sens_Net_heat_flow_ocean_between_surface_and_deep_per_K_of_difference_ZJ_yr_K => 10.0,
        Sens_NF_Avg_life_biomass_yr => 60.0,
        Sens_NF_Speed_of_regrowth_yr => 3.0,
        Sens_Slope_of_effect_of_temp_shifting_DESERT_to_GRASS => 0.4,
        Sens_Slope_temp_vs_glacial_ice_melting => 1.0,
        Sens_Time_in_trunk => 234.638,
        Sens_Time_to_degrade_Kyoto_Flour_yr => 50.0,
        Sens_Time_to_regrow_NF_after_buning_yr => 30.0,
        Sens_TROP_runoff_time => 2000.0,
        Sens_TROP_Time_to_decompose_undisturbed_dead_biomass_yr => 24.0,
        Sensitivity_of_biomass_new_growth_to_CO2_concentration => 1.0,
        Sensitivity_of_convection_to_temp => 2.5,
        Sensitivity_of_evaporation_to_temp => 0.58,
        Sensitivity_of_high_cloud_coverage_to_temp_base => 50.0,
        Sensitivity_of_high_cloud_coverage_to_temp_sens => 50.0,
        Sensitivity_of_low_cloud_coverage_to_temp => 58.0,
        Sensitivity_of_trop_to_humidity => 5.0,
        Slider_for_annual_removal_of_C_from_atm_after_2020_GtC_y => 0.0,
        Slider_for_H2O_slope_hist => 0.0,
        Slider_for_slope_fut => 0.0,
        Slope_btw_Kyoto_Flour_ppt_and_blocking_multiplier => 0.3,
        Slope_btw_Montreal_gases_ppt_and_blocking_multiplier => 0.3,
        Slope_btw_N2O_ppb_and_blocking_multiplier => 0.1,
        Slope_btw_temp_and_permafrost_melting___freezing_base => 1.0,
        Slope_btw_temp_and_permafrost_melting___freezing_sensitivity => 1.0,
        Slope_Effect_Temp_on_NMPP => 2.0,
        Slope_of_effect_of_temp_on_shifting_NF_to_Tundra => 0.1,
        Slope_of_effect_of_temp_on_shifting_TROP_to_NF => 1.0,
        Slope_of_effect_of_temp_shifting_GRASS_to_DESERT => 5.0,
        Slope_of_effect_of_temp_shifting_GRASS_to_NF => 0.1,
        Slope_of_effect_of_temp_shifting_GRASS_to_TROP => 0.2,
        Slope_of_effect_of_temp_shifting_NF_to_GRASS => 0.01,
        Slope_of_effect_of_temp_shifting_NF_to_TROP => 0.2,
        Slope_of_effect_of_temp_shifting_TROP_to_GRASS => 0.05,
        Slope_of_effect_of_temp_shifting_tundra_to_NF => 0.2,
        Slope_of_efffect_of_acidification_on_NMPP => 5.0,
        Slope_temp_eff_on_fire_incidence => 0.1,
        Slope_temp_vs_antarctic_ice_melting => 1.2,
        Slope_temp_vs_Arctic_ice_melting => 0.65,
        Slope_temp_vs_greenland_ice_melting => 0.1,
        Slope_temp_vs_greenland_ice_that_slid_into_the_ocean_melting => 0.71,
        Solar_sine_forcing_amplitude => 0.1,
        Solar_sine_forcing_lift => 0.05,
        Solar_sine_forcing_offset_yr => -3.5,
        Solar_sine_forcing_period_yr => 11.0,
        Stephan_Boltzmann_constant => 5.67037e-8,
        Stratospheric_scattering_experiment_end_year => 3.0e7,
        Stratospheric_scattering_experiment_reduction_from_2015_in_W_m2 => 3.0,
        Switch_0_normal_model_1_dbl_CO2_2_1pct_incr => float(0),
        Switch_btw_historical_CO2_CH4_emissions_or_constant_1history_0constant => float(1),
        SWITCH_for_NATURE_comm_200115_base_1_cut_all_mm_emi_in_2020_2 => float(1),
        SWITCH_future_slope_base_0_plus_5_1_minus_5_2 => float(0),
        SWITCH_h2o_blocked_table_0_linear_1_poly_2 => float(2),
        SWITCH_h2o_poly_dyn_0_equ_1 => float(1),
        SWITCH_nature_rev_0_base_1_steeper_2_less_steep => float(0),
        Switch_to_choose_CO2_CH4_aerosol_other_GHG_0normal_1zero_from_2010_2constant_from_2010 => float(0),
        Switch_to_choose_input_emission_scenario_for_CO2_CH4_and_oth_GHG => float(1),
        Switch_to_drive_model_with_normal_ESCIMO_data__0__CO2e_from_C_Roads__1__or_CO2e_from_CAT_2__or_user_determined_CO2_max_to_find_temp_tipping_point__3_ => 0.0,
        Switch_to_run_experiment_12a_reduction_in_emissions_0_off_1_on => float(0),
        Switch_to_run_experiment_12b_CCS_0_off_1_on => float(0),
        Switch_to_run_experiment_12c_stopping_TROP_deforestation_0_off_1_on => float(0),
        Switch_to_run_experiment_12e_white_surfaces_0_off_1_on => float(0),
        Switch_to_run_NATURE_experiment_CCS_0_off_1_on_0 => float(0),
        Switch_to_run_POLICY_4_Stopping_logging_in_Northern_forests_0_off_1_on => float(0),
        Temp__ocean__deep_in_1850_C => 4.0,
        Temp_atm_1850 => 274.31,
        Temp_gradient_in_surface_degK => 9.7,
        Temp_surface_1850_K => 286.815,
        TEST_Year_in_which_zero_emissions_are_to_be_reached_yr_Remember_to_set_switch_to_9Linear => 2050.0,
        Thickness_of_deep_water_box_1km_to_bottom => 2800.0,
        Thickness_of_intermediate_water_box_800m => 800.0,
        Thickness_of_surface_water_box_100m => 100.0,
        Time_at_which_human_deforestation_is_stopped => 3000.0,
        Time_for_volcanic_aerosols_to_remain_in_the_stratosphere => 1.0,
        Time_in_cold => 6.51772,
        Time_in_deep => 739.89,
        Time_in_intermediate_yr => 211.397,
        Time_in_warm => 26.227,
        Time_to_degrade_Montreal_gases_yr => 30.0,
        Time_to_degrade_N2O_in_atmopshere_yr => 95.0,
        Time_to_deposit_C_in_sediment => 20000.0,
        Time_to_let_shells_form_and_sink_to_sediment_yr => 25.0,
        Time_to_melt_Arctic_ice_at_the_reference_delta_temp => 500.0,
        Time_to_melt_greenland_ice_at_the_reference_delta_temp => 4000.0,
        Time_to_melt_greenland_ice_that_slid_into_the_ocean_at_the_reference_delta_temp => 20.0,
        Time_to_melt_or_freeze_antarctic_ice_at_the_reference_delta_temp => 18000.0,
        Time_to_melt_or_freeze_glacial_ice_at_the_reference_delta_temp => 500.0,
        Time_to_propagate_temperature_change_through_the_volume_of_permafrost_yr => 5.0,
        Time_to_reach_C_equilibrium_between_atmosphere_and_ocean => 18.0,
        Time_to_regrow_GRASS_after_buning_yr => 10.0,
        Time_to_regrow_GRASS_after_deforesting_yr => 80.0,
        Time_to_regrow_NF_after_deforesting_yr => 80.0,
        Time_to_regrow_TROP_after_buning_yr => 30.0,
        Time_to_regrow_TUNDRA_after_buning_yr => 10.0,
        Time_to_regrow_TUNDRA_after_deforesting_yr => 80.0,
        Time_to_smooth_out_temperature_diff_relevant_for_melting_or_freezing_from_1850_yr => 3.0,
        Tipping_point_search_amount_at_peak => 0.0,
        Tipping_point_year_of_end => 210000.0,
        Tipping_point_year_of_start => 500000.0,
        TROP_area_burned_in_1850_Mkm2 => 1.7,
        TROP_area_deforested_in_1850_Mkm2 => 1.0,
        TROP_area_harvested_in_1850_Mkm2 => 0.3,
        TROP_Avg_life_biomass_yr => 60.0,
        TROP_Avg_life_of_building_yr => 20.0,
        TROP_Biomass_locked_in_construction_material_in_1850_GtBiomass => 30.0,
        TROP_clear_cut_fraction => 0.5,
        TROP_Dead_biomass__litter_and_soil_organic_matter_SOM_in_1850_GtBiomass => 160.0,
        TROP_DeadB_and_SOM_densitiy_in_1850_tBiomass_pr_km2 => 8500.0,
        TROP_Fraction_of_construction_waste_burned_0_1 => 0.5,
        TROP_fraction_of_DeadB_and_SOM_being_destroyed_by_clear_cutting => 0.5,
        TROP_fraction_of_DeadB_and_SOM_being_destroyed_by_deforesting => 1.0,
        TROP_fraction_of_DeadB_and_SOM_being_destroyed_by_energy_harvesting => 0.1,
        TROP_fraction_of_DeadB_and_SOM_destroyed_by_natural_fires => 0.0,
        TROP_living_biomass_densitiy_in_1850_tBiomass_pr_km2 => 16500.0,
        TROP_Living_biomass_in_1850_GtBiomass => 370.0,
        TROP_Normal_fire_incidence_fraction_yr => 0.3,
        TROP_Ref_historical_deforestation_pct_yr => 1.0,
        TROP_Slope_temp_eff_on_potential_biomass_per_km2 => -0.5,
        TROP_Speed_of_regrowth_yr => 3.0,
        TUNDRA_area_burned_in_1850_Mkm2 => 2.0,
        TUNDRA_area_deforested_in_1850_Mkm2 => 0.0,
        TUNDRA_area_harvested_in_1850_Mkm2 => 2.5,
        TUNDRA_Avg_life_biomass_yr => 100.0,
        TUNDRA_Avg_life_of_building_yr => 10.0,
        TUNDRA_Biomass_locked_in_construction_material_in_1850_GtBiomass => 1.5,
        TUNDRA_Dead_biomass__litter_and_soil_organic_matter_SOM_in_1850_GtBiomass => 1200.0,
        TUNDRA_DeadB_and_SOM_densitiy_in_1850_tBiomass_pr_km2 => 65000.0,
        TUNDRA_Fraction_of_construction_waste_burned_0_1 => 0.5,
        TUNDRA_fraction_of_DeadB_and_SOM_being_destroyed_by_deforesting => 1.0,
        TUNDRA_fraction_of_DeadB_and_SOM_being_destroyed_by_energy_harvesting => 0.1,
        TUNDRA_fraction_of_DeadB_and_SOM_destroyed_by_natural_fires => 0.0,
        TUNDRA_living_biomass_densitiy_in_1850_tBiomass_pr_km2 => 14500.0,
        TUNDRA_Living_biomass_in_1850_GtBiomass => 300.0,
        TUNDRA_Normal_fire_incidence_fraction_yr => 1.0,
        TUNDRA_Ref_historical_deforestation_pct_yr => 0.0,
        TUNDRA_runoff_time => 2000.0,
        TUNDRA_Speed_of_regrowth_yr => 3.0,
        TUNDRA_Time_to_decompose_undisturbed_dead_biomass_yr => 1000.0,
        UNIT_conversion_1_km3 => 1.0,
        UNIT_conversion_1_yr => 1.0,
        UNIT_conversion_C_to_pH => 1.0,
        UNIT_Conversion_from__km3__km_yr___to_Mkm2_yr => 1.0e-6,
        UNIT_conversion_from_km_to_m => 1000.0,
        UNIT_Conversion_from_km3_to_km2 => 1.0,
        UNIT_Conversion_from_N2O_amount_to_concentration_ppb_MtN2O => 0.305,
        UNIT_conversion_Gm3_to_km3 => 1.0,
        UNIT_conversion_Gt_to_kt => 1.0e6,
        UNIT_conversion_Gt_to_Mt => 1000.0,
        UNIT_conversion_GtBiomass_yr_to_Mkm2_yr => 1000.0,
        UNIT_conversion_GtC_to_MtC => 1000.0,
        UNIT_conversion_GtIce_to_ZJ_melting => 1.0,
        UNIT_conversion_km2___km_to_km3 => 1.0,
        UNIT_conversion_km2_to_Mkm2 => 1.0e6,
        UNIT_conversion_km3_to_Gm3 => 1.0,
        UNIT_conversion_km3_km_to_km2 => 1.0,
        UNIT_conversion_m2_to_km2 => 1.0e6,
        UNIT_conversion_m2_to_Mkm2 => 1.0e12,
        UNIT_conversion_Sv_to_Gm3_yr => 31536.0,
        UNIT_conversion_to_Gm3 => 1.0,
        UNIT_conversion_to_km2_yr => 1.0,
        UNIT_conversion_to_yr => 1.0,
        UNIT_conversion_W_to_ZJ_s => 1.0,
        UNIT_conversion_ymoles___litre_to_dless => 1.0,
        UNIT_conversion_yr_to_dless => 1.0,
        Urban_area_fraction_2000 => 0.004,
        Use_of_GRASS_biomass_for_construction_in_1850_pct => 0.05,
        Use_of_GRASS_biomass_for_energy_in_1850_pct => 1.0,
        Use_of_NF_biomass_for_construction_in_1850_pct => 0.58,
        Use_of_NF_biomass_for_energy_in_1850_pct => 1.09,
        Use_of_TROP_biomass_for_construction_in_1850_pct => 0.48,
        Use_of_TROP_biomass_for_energy_in_1850_pct => 0.07,
        Use_of_TUNDRA_biomass_for_construction_in_1850_pct => 0.05,
        Use_of_TUNDRA_biomass_for_energy_in_1850_pct => 1.0,
        VAES_puls_repetition => 40.0,
        VAES_pulse_duration => 10.0,
        VAES_pulse_height => 1.0,
        Value_of_anthropogenic_aerosol_emissions_during_2015 => 0.225,
        Water_content_of_evaporation_g_kg_per_ZJ_yr => 0.00125,
        Wetlands_area_1850 => 1.0e7,
        When_first_destroyed_yr => float(2020),
        When_methane_hydrates_first_released_yr => float(2020),
        When_to_sample_for_CO2_experiment_yr => float(20000000),
        Yr_to_cut_mm_emi_abrubtly_to_zero_y => 2020.0,
        Zero_C_on_K_scale_K => 273.15,
        Zetta => 1.0e21,
        CO2_concentration_in_1750_ppm => 2.0,
        N2O_ie_N_1750_ppb => 2.0,
        CH4_ie_M_1750_ppb => 2.0,
        LW_Clear_sky_emissions_from_atm_W_m2_in_1850 => 2.0,
        SW_surface_absorption_W_m2_in_1850 => 2.0,
        SW_surface_reflection_W_m2_in_1850 => 2.0,
        C_in_TUNDRA_DeadB_and_soil_in_1850_GtC => 2.0,
        C_in_TUNDRA_LB_in_1850_GtC => 2.0,
        Ga__BB_radiation_less_TOA_radiation_W_m2_in_1850 => 2.0,
        Biomass_new_growing_1850_GtBiomass___yr => 2.0,
        Switch_to_choose_CO2_CH4_aerosol_other_GHG_0normal_1zero_from_202constant_from_2010 => 2.0,
        LW_TOA_radiation_from_atm_to_space_in_1850_W_m2 => 2.0,
      )
      startEquationComponents = []
      begin
        startEquationConstructors = Function[]
        begin
          function generateStartEquations0()
            [
              ifCond1 => false,
              ifCond2 => false,
              Antarctic_ice_volume_km3 => 3.0e7,
              Arctic_ice__on_sea__area_km2 => 1.34e7,
              C_in_atmosphere_GtC => 600.0,
              C_in_atmosphere_in_form_of_CH4 => 1.69,
              C_in_cold_surface_water_GtC => Carbon_in_cold_ocean_0_to_100m_1850_GtC,
              C_in_cold_water_trunk_downwelling_GtC => Carbon_in_cold_ocean_trunk_100m_to_bottom_1850_GtC,
              C_in_deep_water_volume_1km_to_bottom_GtC => Carbon_in_ocean_deep_1k_to_bottom_ocean_1850_GtC,
              C_in_intermediate_upwelling_water_100m_to_1km_GtC => Carbon_in_ocean_upwelling_100m_to_1km_1850_GtC,
              C_in_permafrost_in_form_of_CH4 => 1200.0,
              C_in_sediment => 3.0e9,
              C_in_warm_surface_water_GtC => Carbon_in_warm_ocean_0_to_100m_1850_GtC,
              Cold_surface_water_volume_Gm3 => Volume_cold_ocean_0_to_100m,
              Cold_water_volume_downwelling_Gm3 => Volume_cold_ocean_downwelling_100m_to_bottom,
              Cumulative_antarctic_ice_volume_loss_GtIce => 0.0,
              Cumulative_C_released_from_permafrost_as_either_CH4_or_CO2_in_GtC_yr => 0.0,
              Cumulative_carbon_captured_and_stored_GtC => 0.0,
              Cumulative_carbon_removed_from_atm_for_nature_May_2020 => 0.0,
              Cumulative_flow_of_C_to_biomass_since_1850_GtC => 0.0,
              Cumulative_glacial_ice_volume_loss_GtIce => 0.0,
              Cumulative_Greenland_ice_volume_loss_GtIce => 0.0,
              Cumulative_heat_to_atm_ZJ => 0.0,
              Cumulative_ocean_volume_increase_due_to_ice_melting_km3 => 0.0,
              Cumulative_release_of_C_from_permafrost_GtC => 0.0,
              Deep_water_volume_1km_to_4km_Gm3 => Volume_ocean_deep_1km_to_bottom,
              DESERT_Mkm2 => 25.4,
              Fossil_fuel_reserves_in_ground_GtC => 6000.0,
              Glacial_ice_volume_km3 => 167000.0,
              GRASS_area_burnt_Mkm2 => 1.0,
              GRASS_area_harvested_Mkm2 => 2.5,
              GRASS_Biomass_locked_in_construction_material_GtBiomass => 1.5,
              GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass => 1200.0,
              GRASS_deforested_Mkm2 => 0.5,
              GRASS_Living_biomass_GtBiomass => 310.0,
              GRASS_potential_area_Mkm2 => 22.5,
              Greenland_ice_volume_on_Greenland_km3 => 2.93e6,
              Greenland_ice_volume_that_slid_into_the_ocean_km3 => 0.0,
              Heat_in_atmosphere_ZJ => 1025.67,
              Heat_in_deep_ZJ => 1.9532e6,
              Heat_in_surface => 25000.0,
              Intermediate_upwelling_water_volume_100m_to_1km_Gm3 => Volume_ocean_upwelling_100m_to_1km,
              Kyoto_Flour_gases_in_atm => 0.0,
              Montreal_gases_in_atm => 0.0,
              N2O_in_atmosphere_MtN2O => 900.0,
              NATURE_Cumulative_CCS_GtC => 0.0,
              NF_area_burnt_Mkm2 => 2.5,
              NF_area_clear_cut_Mkm2 => 1.0,
              NF_area_deforested_Mkm2 => 0.0,
              NF_area_harvested_Mkm2 => 1.0,
            ]
          end
          push!(startEquationConstructors, generateStartEquations0)
        end
        begin
          function generateStartEquations1()
            [
              NF_Biomass_locked_in_construction_material_GtBiomass => 3.0,
              NF_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass => 330.0,
              NF_Living_biomass_GtBiomass => 115.0,
              NF_potential_area_Mkm2 => 17.0,
              Sum_C_absorbed_by_ocean_GtC => 0.0,
              Sum_heat_to_deep_ocean => 0.0,
              Sum_heat_to_deep_ocean_btw_72_and_08 => 0.0,
              Sum_heat_to_surface_ocean_btw_72_and_08 => 0.0,
              Sum_heat_to_surface_ocean_ZJ => 0.0,
              Sum_man_made_CO2_emissions_GtC => 0.0,
              Sum_net_C_to_atm => 0.0,
              TROP_area_burnt_Mkm2 => 1.7,
              TROP_area_clear_cut_Mkm2 => 0.3,
              TROP_area_deforested_Mkm2 => 1.0,
              TROP_area_harvested_Mkm2 => 0.3,
              TROP_Biomass_locked_in_construction_material_GtBiomass => 30.0,
              TROP_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass => 160.0,
              TROP_Living_biomass_GtBiomass => 370.0,
              TROP_potential_area_Mkm2 => 25.0,
              TUNDRA_area_burnt_Mkm2 => 2.0,
              TUNDRA_area_harvested_Mkm2 => 2.5,
              TUNDRA_Biomass_locked_in_construction_material_GtBiomass => 1.5,
              TUNDRA_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass => 1200.0,
              TUNDRA_deforested_Mkm2 => 0.0,
              TUNDRA_Living_biomass_GtBiomass => 300.0,
              Tundra_potential_area_Mkm2 => 22.5,
              Volcanic_aerosols_in_stratosphere => 0.0,
              Warm_surface_water_volume_Gm3 => Volume_warm_ocean_0_to_100m,
              Wetlands_area => 1.0e7,
              Aerosol_anthropogenic_emissions_in_2010 => 0.0,
              CO2_emissions_in_2010 => 0.0,
              CO2_ppm_value_at_When_to_sample => MODEL_CO2_concentration_in_atmosphere2_ppm,
              CO4_emissions_in_2010 => 0.0,
              Greenland_slide_experiment_end_condition => 0.0,
              Kyoto_Flour_concentration_in_1970_ppt => 0.0,
              Kyoto_Flour_emissions_RCPs_JR_in_2010 => 0.0,
              Montreal_gases_concentration_in_1970_ppt => 0.0,
              Montreal_gases_emissions_RCPs_JR_in_2010 => 0.0,
              N20_emissions_RCPs_JR_in_2010 => 0.0,
              Tipping_point_search_amount_at_start => 12.0,
              Arctic_land_surface_temp_anomaly_compared_to_1850 => Temp_surface_anomaly_compared_to_1850_degC,
              Biological_removal_of_C_from_WSW_GtC_per_yr => Net_marine_primary_production_NMPP_GtC_pr_yr,
              Effect_of_temp_on_permafrost_melting_dmnl => 1.0 + Slope_btw_temp_and_permafrost_melting___freezing * (Temp_diff_relevant_for_melting_or_freezing_from_1850 / 4.0 - 1.0),
              Temp_diff_relevant_for_melting_or_freezing_arctic_ice_from_1850 => Temp_surface_anomaly_compared_to_1850_degC,
              Temp_diff_relevant_for_melting_or_freezing_from_1850 => Temp_surface_C - 13.66500000000002,
              yr_on_yr_change_in_C_in_atm_GtC_yr => C_in_atmosphere_GtC - C_in_atm_1_yr_ago_GtC,
              C_in_ocean_1_yr_ago_GtC => Total_carbon_in_ocean_GtC,
              C_in_ocean_1_yr_ago_GtC_LV1 => Total_carbon_in_ocean_GtC,
              C_in_ocean_1_yr_ago_GtC_LV2 => Total_carbon_in_ocean_GtC,
              Atmos_heat_used_for_melting_last_year_1_yr_LV => 0.0,
            ]
          end
          push!(startEquationConstructors, generateStartEquations1)
        end
        begin
          function generateStartEquations2()
            [
              Ocean_heat_used_for_melting_last_year_ZJ_yr_LV => 0.0,
              C_in_atm_1_yr_ago_GtC_LV3 => C_in_atm_1_yr_ago_GtC_DL * C_in_atmosphere_GtC,
              C_in_atm_1_yr_ago_GtC_LV2 => C_in_atm_1_yr_ago_GtC_LV3,
              C_in_atm_1_yr_ago_GtC_LV1 => C_in_atm_1_yr_ago_GtC_LV3,
              All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_LV3 => All_C_taken_out_due_to_change_in_land_use_GtC * All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_DL,
              All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_LV2 => All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_LV3,
              All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_LV1 => All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_LV3,
              Model_N2O_concentration_in_1850_ppb => 0.0,
              CO2_concentration_in_1850_ppm => 0.0,
              Incoming_solar_in_1850_ZJ_yr => 0.0,
              C_in_atmosphere_GtC_in_1850 => 0.0,
              C_in_biomass_in_1850_GtC => 0.0,
              Total_carbon_in_ocean_GtC_in_1850 => 0.0,
              Temp_ocean_deep_1850_degC => 0.0,
              init_ph_in_cold_water => 0.0,
              Humidity_of_atmosphere_in_1850_g_kg => 0.0,
              LW_TOA_radiation_from_atm_to_space_in_1850 => 0.0,
              Temp__ocean__surface_in_1850_C => 0.0,
              Fraction_blocked_by_ALL_GHG_in_1850 => 0.0,
              Fraction_blocked_CO2_in_1850 => 0.0,
              Fraction_blocked_CH4_in_1850 => 0.0,
              Fraction_blocked_othGHG_in_1850 => 0.0,
              init_C_in_GRASS => 0.0,
              init_C_in_NF => 0.0,
              init_C_in_TROP => 0.0,
              init_C_in_TUNDRA => 0.0,
              Fossil_fuel_reserves_in_ground_1850_GtC => 0.0,
              Time => 0.0,
              Aerosol_anthropogenic_emissions_in_2010 => 0.0,
              CO2_emissions_in_2010 => 0.0,
              CO2_ppm_value_at_When_to_sample => 0.0,
              CO4_emissions_in_2010 => 0.0,
              Greenland_slide_experiment_end_condition => 0.0,
              Kyoto_Flour_concentration_in_1970_ppt => 0.0,
              Kyoto_Flour_emissions_RCPs_JR_in_2010 => 0.0,
              Montreal_gases_concentration_in_1970_ppt => 0.0,
              Montreal_gases_emissions_RCPs_JR_in_2010 => 0.0,
              N20_emissions_RCPs_JR_in_2010 => 0.0,
              Tipping_point_search_amount_at_start => 0.0,
              combi_E3_SC_1_CO2_GtC_yr_u => 0.0,
              var"combi_E3_SC_1_CO2_GtC_yr_y[1]" => 0.0,
              E3_SC_1_CO2_GtC_yr => 0.0,
              combi_E3_SC_1_CH4_GtC_yr_u => 0.0,
              var"combi_E3_SC_1_CH4_GtC_yr_y[1]" => 0.0,
              E3_SC_1_CH4_GtC_yr => 0.0,
              combi_E3_SC_1_N2O_Mt_yr_u => 0.0,
              var"combi_E3_SC_1_N2O_Mt_yr_y[1]" => 0.0,
              E3_SC_1_N2O_Mt_yr => 0.0,
              combi_E3_SC_1_Kyoto_F_kt_yr_u => 0.0,
              var"combi_E3_SC_1_Kyoto_F_kt_yr_y[1]" => 0.0,
            ]
          end
          push!(startEquationConstructors, generateStartEquations2)
        end
        begin
          function generateStartEquations3()
            [
              E3_SC_1_Kyoto_F_kt_yr => 0.0,
              combi_E3_SC_1_Montreal_gases_kt_yr_u => 0.0,
              var"combi_E3_SC_1_Montreal_gases_kt_yr_y[1]" => 0.0,
              E3_SC_1_Montreal_gases_kt_yr => 0.0,
              combi_E3_SC_2_CO2_GtC_yr_u => 0.0,
              var"combi_E3_SC_2_CO2_GtC_yr_y[1]" => 0.0,
              E3_SC_2_CO2_GtC_yr => 0.0,
              combi_E3_SC_2_CH4_GtC_yr_u => 0.0,
              var"combi_E3_SC_2_CH4_GtC_yr_y[1]" => 0.0,
              E3_SC_2_CH4_GtC_yr => 0.0,
              combi_E3_SC_2_N2O_Mt_yr_u => 0.0,
              var"combi_E3_SC_2_N2O_Mt_yr_y[1]" => 0.0,
              E3_SC_2_N2O_Mt_yr => 0.0,
              combi_E3_SC_2_Kyoto_F_kt_yr_u => 0.0,
              var"combi_E3_SC_2_Kyoto_F_kt_yr_y[1]" => 0.0,
              E3_SC_2_Kyoto_F_kt_yr => 0.0,
              combi_E3_SC_2_Montreal_gases_kt_yr_u => 0.0,
              var"combi_E3_SC_2_Montreal_gases_kt_yr_y[1]" => 0.0,
              E3_SC_2_Montreal_gases_kt_yr => 0.0,
              combi_E3_SC_3_CO2_GtC_yr_u => 0.0,
              var"combi_E3_SC_3_CO2_GtC_yr_y[1]" => 0.0,
              E3_SC_3_CO2_GtC_yr => 0.0,
              combi_E3_SC_3_CH4_GtC_yr_u => 0.0,
              var"combi_E3_SC_3_CH4_GtC_yr_y[1]" => 0.0,
              E3_SC_3_CH4_GtC_yr => 0.0,
              combi_E3_SC_3_N2O_Mt_yr_u => 0.0,
              var"combi_E3_SC_3_N2O_Mt_yr_y[1]" => 0.0,
              E3_SC_3_N2O_Mt_yr => 0.0,
              combi_E3_SC_3_Kyoto_F_kt_yr_u => 0.0,
              var"combi_E3_SC_3_Kyoto_F_kt_yr_y[1]" => 0.0,
              E3_SC_3_Kyoto_F_kt_yr => 0.0,
              combi_E3_SC_3_Montreal_gases_kt_yr_u => 0.0,
              var"combi_E3_SC_3_Montreal_gases_kt_yr_y[1]" => 0.0,
              E3_SC_3_Montreal_gases_kt_yr => 0.0,
              combi_E3_SC_4_CO2_GtC_yr_u => 0.0,
              var"combi_E3_SC_4_CO2_GtC_yr_y[1]" => 0.0,
              E3_SC_4_CO2_GtC_yr => 0.0,
              combi_E3_SC_4_CH4_GtC_yr_u => 0.0,
              var"combi_E3_SC_4_CH4_GtC_yr_y[1]" => 0.0,
              E3_SC_4_CH4_GtC_yr => 0.0,
              combi_E3_SC_4_N2O_Mt_yr_u => 0.0,
              var"combi_E3_SC_4_N2O_Mt_yr_y[1]" => 0.0,
              E3_SC_4_N2O_Mt_yr => 0.0,
              combi_E3_SC_4_Kyoto_F_kt_yr_u => 0.0,
              var"combi_E3_SC_4_Kyoto_F_kt_yr_y[1]" => 0.0,
              E3_SC_4_Kyoto_F_kt_yr => 0.0,
              combi_E3_SC_4_Montreal_gases_kt_yr_u => 0.0,
              var"combi_E3_SC_4_Montreal_gases_kt_yr_y[1]" => 0.0,
              E3_SC_4_Montreal_gases_kt_yr => 0.0,
              combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_MtCH4_yr_u => 0.0,
            ]
          end
          push!(startEquationConstructors, generateStartEquations3)
        end
        begin
          function generateStartEquations4()
            [
              var"combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_MtCH4_yr_y[1]" => 0.0,
              Emissions_of_anthro_CH4_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_MtCH4_yr => 0.0,
              combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP3_MtCH4_yr_u => 0.0,
              var"combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP3_MtCH4_yr_y[1]" => 0.0,
              Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP3_MtCH4_yr => 0.0,
              combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP45_MtCH4_yr_u => 0.0,
              var"combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP45_MtCH4_yr_y[1]" => 0.0,
              Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP45_MtCH4_yr => 0.0,
              combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP6_MtCH4_yr_u => 0.0,
              var"combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP6_MtCH4_yr_y[1]" => 0.0,
              Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP6_MtCH4_yr => 0.0,
              combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP85_MtCH4_yr_u => 0.0,
              var"combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP85_MtCH4_yr_y[1]" => 0.0,
              Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP85_MtCH4_yr => 0.0,
              combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_GtC_yr_u => 0.0,
              var"combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_GtC_yr_y[1]" => 0.0,
              Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_GtC_yr => 0.0,
              combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_GtC_yr_u => 0.0,
              var"combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_GtC_yr_y[1]" => 0.0,
              Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_GtC_yr => 0.0,
              combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_hist_with_JR__large_scale_CCS__exp_GtC_yr_u => 0.0,
              var"combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_hist_with_JR__large_scale_CCS__exp_GtC_yr_y[1]" => 0.0,
              Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_hist_with_JR__large_scale_CCS__exp_GtC_yr => 0.0,
              combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP3_GtC_yr_u => 0.0,
              var"combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP3_GtC_yr_y[1]" => 0.0,
              Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP3_GtC_yr => 0.0,
              combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP45_GtC_yr_u => 0.0,
              var"combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP45_GtC_yr_y[1]" => 0.0,
              Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP45_GtC_yr => 0.0,
              combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP6_GtC_yr_u => 0.0,
              var"combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP6_GtC_yr_y[1]" => 0.0,
              Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP6_GtC_yr => 0.0,
              combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP85_GtC_yr_u => 0.0,
              var"combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP85_GtC_yr_y[1]" => 0.0,
              Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP85_GtC_yr => 0.0,
              combi_Emissions_of_anthro_N20_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_MtN2O_yr_u => 0.0,
              var"combi_Emissions_of_anthro_N20_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_MtN2O_yr_y[1]" => 0.0,
              Emissions_of_anthro_N20_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_MtN2O_yr => 0.0,
              combi_Emissions_of_anthro_N20_with_JR_2052_shape_MtN20_yr_u => 0.0,
              var"combi_Emissions_of_anthro_N20_with_JR_2052_shape_MtN20_yr_y[1]" => 0.0,
              Emissions_of_anthro_N20_with_JR_2052_shape_MtN20_yr => 0.0,
              combi_Emissions_of_Kyoto_Flour_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_kt_yr_u => 0.0,
              var"combi_Emissions_of_Kyoto_Flour_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_kt_yr_y[1]" => 0.0,
              Emissions_of_Kyoto_Flour_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_kt_yr => 0.0,
              combi_Emissions_of_Kyoto_Flour_with_JR_2052_shape_kt_yr_u => 0.0,
              var"combi_Emissions_of_Kyoto_Flour_with_JR_2052_shape_kt_yr_y[1]" => 0.0,
              Emissions_of_Kyoto_Flour_with_JR_2052_shape_kt_yr => 0.0,
              combi_Emissions_of_Montreal_gases_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_kt_yr_u => 0.0,
              var"combi_Emissions_of_Montreal_gases_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_kt_yr_y[1]" => 0.0,
              Emissions_of_Montreal_gases_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_kt_yr => 0.0,
            ]
          end
          push!(startEquationConstructors, generateStartEquations4)
        end
        begin
          function generateStartEquations5()
            [
              combi_Emissions_of_Montreal_gases_with_JR_2052_shape_kt_yr_u => 0.0,
              var"combi_Emissions_of_Montreal_gases_with_JR_2052_shape_kt_yr_y[1]" => 0.0,
              Emissions_of_Montreal_gases_with_JR_2052_shape_kt_yr => 0.0,
              combi_CH4_emissions_from_CO2e_C_Roads_u => 0.0,
              var"combi_CH4_emissions_from_CO2e_C_Roads_y[1]" => 0.0,
              CH4_emissions_from_CO2e_C_Roads => 0.0,
              combi_CH4_emissions_from_CO2e_CAT_u => 0.0,
              var"combi_CH4_emissions_from_CO2e_CAT_y[1]" => 0.0,
              CH4_emissions_from_CO2e_CAT => 0.0,
              combi_CH4_emissions_pct_contribution_to_Total_CO2e_u => 0.0,
              var"combi_CH4_emissions_pct_contribution_to_Total_CO2e_y[1]" => 0.0,
              CH4_emissions_pct_contribution_to_Total_CO2e => 0.0,
              combi_CO2_emissions_from_CO2e_C_Roads_u => 0.0,
              var"combi_CO2_emissions_from_CO2e_C_Roads_y[1]" => 0.0,
              CO2_emissions_from_CO2e_C_Roads => 0.0,
              combi_CO2_emissions_from_CO2e_CAT_u => 0.0,
              var"combi_CO2_emissions_from_CO2e_CAT_y[1]" => 0.0,
              CO2_emissions_from_CO2e_CAT => 0.0,
              combi_CO2_emissions_pct_contribution_to_Total_CO2e_u => 0.0,
              var"combi_CO2_emissions_pct_contribution_to_Total_CO2e_y[1]" => 0.0,
              CO2_emissions_pct_contribution_to_Total_CO2e => 0.0,
              combi_Historical_aerosol_emissions_anthro_u => 0.0,
              var"combi_Historical_aerosol_emissions_anthro_y[1]" => 0.0,
              Historical_aerosol_emissions_anthro => 0.0,
              combi_Historical_forcing_from_solar_insolation_W_m2_u => 0.0,
              var"combi_Historical_forcing_from_solar_insolation_W_m2_y[1]" => 0.0,
              Historical_forcing_from_solar_insolation_W_m2 => 0.0,
              combi_Historical_aerosol_forcing_volcanic_u => 0.0,
              var"combi_Historical_aerosol_forcing_volcanic_y[1]" => 0.0,
              Historical_aerosol_forcing_volcanic => 0.0,
              combi_OGHG_Kyoto_Flour_emi_rcp3_u => 0.0,
              var"combi_OGHG_Kyoto_Flour_emi_rcp3_y[1]" => 0.0,
              OGHG_Kyoto_Flour_emi_rcp3 => 0.0,
              combi_OGHG_Kyoto_Flour_emi_rcp45_u => 0.0,
              var"combi_OGHG_Kyoto_Flour_emi_rcp45_y[1]" => 0.0,
              OGHG_Kyoto_Flour_emi_rcp45 => 0.0,
              combi_OGHG_Kyoto_Flour_emi_rcp6_u => 0.0,
              var"combi_OGHG_Kyoto_Flour_emi_rcp6_y[1]" => 0.0,
              OGHG_Kyoto_Flour_emi_rcp6 => 0.0,
              combi_OGHG_Kyoto_Flour_emi_rcp85_u => 0.0,
              var"combi_OGHG_Kyoto_Flour_emi_rcp85_y[1]" => 0.0,
              OGHG_Kyoto_Flour_emi_rcp85 => 0.0,
              combi_Kyoto_Flour_emissions_from_CO2e_C_Roads_u => 0.0,
              var"combi_Kyoto_Flour_emissions_from_CO2e_C_Roads_y[1]" => 0.0,
              Kyoto_Flour_emissions_from_CO2e_C_Roads => 0.0,
              combi_Kyoto_Flour_emissions_from_CO2e_CAT_u => 0.0,
              var"combi_Kyoto_Flour_emissions_from_CO2e_CAT_y[1]" => 0.0,
              Kyoto_Flour_emissions_from_CO2e_CAT => 0.0,
              combi_Kyoto_Flour_emissions_pct_contribution_to_Total_CO2e_u => 0.0,
              var"combi_Kyoto_Flour_emissions_pct_contribution_to_Total_CO2e_y[1]" => 0.0,
            ]
          end
          push!(startEquationConstructors, generateStartEquations5)
        end
        begin
          function generateStartEquations6()
            [
              Kyoto_Flour_emissions_pct_contribution_to_Total_CO2e => 0.0,
              combi_OGHG_Montreal_gases_emi_rcp3_u => 0.0,
              var"combi_OGHG_Montreal_gases_emi_rcp3_y[1]" => 0.0,
              OGHG_Montreal_gases_emi_rcp3 => 0.0,
              combi_OGHG_Montreal_gases_emi_rcp45_u => 0.0,
              var"combi_OGHG_Montreal_gases_emi_rcp45_y[1]" => 0.0,
              OGHG_Montreal_gases_emi_rcp45 => 0.0,
              combi_OGHG_Montreal_gases_emi_rcp6_u => 0.0,
              var"combi_OGHG_Montreal_gases_emi_rcp6_y[1]" => 0.0,
              OGHG_Montreal_gases_emi_rcp6 => 0.0,
              combi_OGHG_Montreal_gases_emi_rcp85_u => 0.0,
              var"combi_OGHG_Montreal_gases_emi_rcp85_y[1]" => 0.0,
              OGHG_Montreal_gases_emi_rcp85 => 0.0,
              combi_othGHG_N20_man_made_emissions_rcp3_u => 0.0,
              var"combi_othGHG_N20_man_made_emissions_rcp3_y[1]" => 0.0,
              othGHG_N20_man_made_emissions_rcp3 => 0.0,
              combi_othGHG_N20_man_made_emissions_rcp45_u => 0.0,
              var"combi_othGHG_N20_man_made_emissions_rcp45_y[1]" => 0.0,
              othGHG_N20_man_made_emissions_rcp45 => 0.0,
              combi_othGHG_N20_man_made_emissions_rcp6_u => 0.0,
              var"combi_othGHG_N20_man_made_emissions_rcp6_y[1]" => 0.0,
              othGHG_N20_man_made_emissions_rcp6 => 0.0,
              combi_othGHG_N20_man_made_emissions_rcp85_u => 0.0,
              var"combi_othGHG_N20_man_made_emissions_rcp85_y[1]" => 0.0,
              othGHG_N20_man_made_emissions_rcp85 => 0.0,
              combi_RCP_3_CO2_concentration_1850_2100_ppm_u => 0.0,
              var"combi_RCP_3_CO2_concentration_1850_2100_ppm_y[1]" => 0.0,
              RCP_3_CO2_concentration_1850_2100_ppm => 0.0,
              combi_RCP_45_CO2_concentration_1850_2100_ppm_u => 0.0,
              var"combi_RCP_45_CO2_concentration_1850_2100_ppm_y[1]" => 0.0,
              RCP_45_CO2_concentration_1850_2100_ppm => 0.0,
              combi_RCP_6_CO2_concentration_1850_2100_ppm_u => 0.0,
              var"combi_RCP_6_CO2_concentration_1850_2100_ppm_y[1]" => 0.0,
              RCP_6_CO2_concentration_1850_2100_ppm => 0.0,
              combi_RCP_85_CO2_concentration_1850_2100_ppm_u => 0.0,
              var"combi_RCP_85_CO2_concentration_1850_2100_ppm_y[1]" => 0.0,
              RCP_85_CO2_concentration_1850_2100_ppm => 0.0,
              combi_Montreal_gases_emissions_from_CO2e_C_Roads_u => 0.0,
              var"combi_Montreal_gases_emissions_from_CO2e_C_Roads_y[1]" => 0.0,
              Montreal_gases_emissions_from_CO2e_C_Roads => 0.0,
              combi_Montreal_gases_emissions_from_CO2e_CAT_u => 0.0,
              var"combi_Montreal_gases_emissions_from_CO2e_CAT_y[1]" => 0.0,
              Montreal_gases_emissions_from_CO2e_CAT => 0.0,
              combi_Montreal_gases_emissions_pct_contribution_to_Total_CO2e_u => 0.0,
              var"combi_Montreal_gases_emissions_pct_contribution_to_Total_CO2e_y[1]" => 0.0,
              Montreal_gases_emissions_pct_contribution_to_Total_CO2e => 0.0,
              combi_N2O_man_made_emissions_from_CO2e_C_Roads_u => 0.0,
              var"combi_N2O_man_made_emissions_from_CO2e_C_Roads_y[1]" => 0.0,
              N2O_man_made_emissions_from_CO2e_C_Roads => 0.0,
              combi_N2O_man_made_emissions_from_CO2e_CAT_u => 0.0,
            ]
          end
          push!(startEquationConstructors, generateStartEquations6)
        end
        begin
          function generateStartEquations7()
            [
              var"combi_N2O_man_made_emissions_from_CO2e_CAT_y[1]" => 0.0,
              N2O_man_made_emissions_from_CO2e_CAT => 0.0,
              combi_N2O_emissions_pct_contribution_to_Total_CO2e_u => 0.0,
              var"combi_N2O_emissions_pct_contribution_to_Total_CO2e_y[1]" => 0.0,
              N2O_emissions_pct_contribution_to_Total_CO2e => 0.0,
              combi_Sea_level_rise_history_mm_u => 0.0,
              var"combi_Sea_level_rise_history_mm_y[1]" => 0.0,
              Sea_level_rise_history_mm => 0.0,
              combi_All_Human_activity_emissions_GtCO2e_yr_Base_for_tipping_point_search_u => 0.0,
              var"combi_All_Human_activity_emissions_GtCO2e_yr_Base_for_tipping_point_search_y[1]" => 0.0,
              combi_Arctic_freezing_cutoff_u => 0.0,
              var"combi_Arctic_freezing_cutoff_y[1]" => 0.0,
              combi_Blocked_by_H20_hist_Table_lookup_u => 0.0,
              var"combi_Blocked_by_H20_hist_Table_lookup_y[1]" => 0.0,
              combi_Blocked_by_H20_Table_lookup_u => 0.0,
              var"combi_Blocked_by_H20_Table_lookup_y[1]" => 0.0,
              combi_Effect_of_heat_in_atm_on_melting_ice__cut_off__u => 0.0,
              var"combi_Effect_of_heat_in_atm_on_melting_ice__cut_off__y[1]" => 0.0,
              combi_Exp_12a_reduction_in_emissions_LOOKUP_u => 0.0,
              var"combi_Exp_12a_reduction_in_emissions_LOOKUP_y[1]" => 0.0,
              combi_EXP_12b_CCS_from_2015_u => 0.0,
              var"combi_EXP_12b_CCS_from_2015_y[1]" => 0.0,
              combi_EXP_12e_white_surfaces_ease_in_u => 0.0,
              var"combi_EXP_12e_white_surfaces_ease_in_y[1]" => 0.0,
              combi_Fraction_blocked_by_CH4_spectrum_u => 0.0,
              var"combi_Fraction_blocked_by_CH4_spectrum_y[1]" => 0.0,
              combi_Fraction_blocked_by_CO2_spectrum_u => 0.0,
              var"combi_Fraction_blocked_by_CO2_spectrum_y[1]" => 0.0,
              combi_Future_shape_of_anthropogenic_aerosol_emissions_u => 0.0,
              var"combi_Future_shape_of_anthropogenic_aerosol_emissions_y[1]" => 0.0,
              combi_Melting_constraint_from_the_heat_in__ocean__surface_reservoir_u => 0.0,
              var"combi_Melting_constraint_from_the_heat_in__ocean__surface_reservoir_y[1]" => 0.0,
              combi_Melting_constraint_from_the_heat_in_atmosphere_reservoir_fraction_u => 0.0,
              var"combi_Melting_constraint_from_the_heat_in_atmosphere_reservoir_fraction_y[1]" => 0.0,
              combi_Melting_restraint_for_permafrost_from_heat_in_atmophere_u => 0.0,
              var"combi_Melting_restraint_for_permafrost_from_heat_in_atmophere_y[1]" => 0.0,
              combi_NATURE_CCS_removal_experiment_multiplier_u => 0.0,
              var"combi_NATURE_CCS_removal_experiment_multiplier_y[1]" => 0.0,
              combi_NF_clear_cut_fraction_u => 0.0,
              var"combi_NF_clear_cut_fraction_y[1]" => 0.0,
              combi_NF_usage_cutoff_u => 0.0,
              var"combi_NF_usage_cutoff_y[1]" => 0.0,
              combi_Permafrost_melting_cutoff_u => 0.0,
              var"combi_Permafrost_melting_cutoff_y[1]" => 0.0,
              combi_RCPFossil_fuel_usage_cutoff_u => 0.0,
              var"combi_RCPFossil_fuel_usage_cutoff_y[1]" => 0.0,
              combi_Snowball_earth_cutoff_u => 0.0,
              var"combi_Snowball_earth_cutoff_y[1]" => 0.0,
              combi_Thermal_expansion_deep_in_1850_pct_u => 0.0,
              var"combi_Thermal_expansion_deep_in_1850_pct_y[1]" => 0.0,
            ]
          end
          push!(startEquationConstructors, generateStartEquations7)
        end
        begin
          function generateStartEquations8()
            [
              combi_Thermal_expansion_deep_pct_u => 0.0,
              var"combi_Thermal_expansion_deep_pct_y[1]" => 0.0,
              combi_Thermal_expansion_surface_in_1850_pct_u => 0.0,
              var"combi_Thermal_expansion_surface_in_1850_pct_y[1]" => 0.0,
              combi_Thermal_expansion_surface_pct_u => 0.0,
              var"combi_Thermal_expansion_surface_pct_y[1]" => 0.0,
              combi_TROP_deforestation_cutoff_u => 0.0,
              var"combi_TROP_deforestation_cutoff_y[1]" => 0.0,
              combi_TROP_deforestation_cutoff_effect_u => 0.0,
              var"combi_TROP_deforestation_cutoff_effect_y[1]" => 0.0,
              combi_TROP_deforestion_multiplier_wrt_2000_u => 0.0,
              var"combi_TROP_deforestion_multiplier_wrt_2000_y[1]" => 0.0,
              combi_Urbanzation_Effect_on_biomass_use_u => 0.0,
              var"combi_Urbanzation_Effect_on_biomass_use_y[1]" => 0.0,
              combi_Population_Lookup_bn_u => 0.0,
              var"combi_Population_Lookup_bn_y[1]" => 0.0,
              aux_1____Temp_gradient_minus_1___slope_ => 0.0,
              Actual_time_to_degrade_all_GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_yr => 0.0,
              Actual_time_to_degrade_all_NF_Dead_biomass__litter_and_soil_organic_matter_SOM_yr => 0.0,
              Actual_time_to_degrade_all_TROP_Dead_biomass__litter_and_soil_organic_matter_SOM_yr => 0.0,
              Actual_time_to_degrade_all_TUNDRA_Dead_biomass__litter_and_soil_organic_matter_SOM_yr => 0.0,
              Aerosol_anthropogenic_emissions => 0.0,
              Albedo_Antartic => 0.0,
              Albedo_glacier => 0.0,
              Albedo_land_biomes => 0.0,
              Albedo_ocean_with_arctic_ice_changes => 0.0,
              Albedo_URBAN => 0.0,
              All_C_taken_out_due_to_change_in_land_use_GtC => 0.0,
              All_CH4_emissions_GtC_yr => 0.0,
              ALL_clouds_net_effect__pos_warming__neg_cooling__W_m2 => 0.0,
              All_Human_activity_emissions_GtCO2e_yr_Base_for_tipping_point_search => 0.0,
              All_Human_activity_emissions_GtCO2e_yr => 0.0,
              All_N2O_emissions_MtN2O_yr => 0.0,
              Annual_flux_of_C_to_biomass_GtC_pr_yr => 0.0,
              Annual_glacial_ice_losing__pos__or_gaining__neg__GtIce_yr => 0.0,
              Annual_release_of_C_from_permafrost_GtC_y => 0.0,
              Antarctic_ice_area_decrease_Mkm2_pr_yr => 0.0,
              Antarctic_ice_area_increase_Mkm2_pr_yr => 0.0,
              Antarctic_ice_area_km2 => 0.0,
              Antarctic_ice_freezing_km3_yr => 0.0,
              Antarctic_ice_losing__pos__or_gaining__neg__GtIce_yr => 0.0,
              Antarctic_ice_melting__pos__or_freezing__neg__km3_yr => 0.0,
              Antarctic_ice_melting_as_water_km3_yr => 0.0,
              Antarctic_ice_melting_km3_yr => 0.0,
              Anthropogenic_aerosol_forcing => 0.0,
              Arctic_as_fraction_of_ocean => 0.0,
              Arctic_freezing_cutoff => 0.0,
              Arctic_ice_area_max_km2 => 0.0,
              Arctic_ice_area_Mkm2 => 0.0,
              Arctic_ice_melting__pos__or_freezing__neg__km2_yr => 0.0,
            ]
          end
          push!(startEquationConstructors, generateStartEquations8)
        end
        begin
          function generateStartEquations9()
            [
              Area_covered_by_high_clouds => 0.0,
              Area_covered_by_low_clouds => 0.0,
              Area_equivalent_of_linear_retreat_km2_yr => 0.0,
              Area_of_earth_Mkm2 => 0.0,
              Area_of_land_Mkm2 => 0.0,
              Atmos_heat_used_for_melting_1_yr => 0.0,
              Avg_C_concentration_in_top_layer => 0.0,
              Avg_CC_in_ocean_top_layer_ymoles_per_litre => 0.0,
              Avg_CO2_conc_in_ocean_top_layer_in_ppm => 0.0,
              Avg_earths_surface_albedo => 0.0,
              Avg_thickness_Antarctic_km => 0.0,
              Avg_thickness_glacier_km => 0.0,
              Avg_volcanic_activity_GtC_yr => 0.0,
              Barren_land_Mkm2 => 0.0,
              BARREN_land_normal_albedo_Mkm2 => 0.0,
              BARREN_land_white_Mkm2 => 0.0,
              BB_radiation_at_atm_temp_in_atm_W_m2 => 0.0,
              BB_radiation_at_surface_temp_ZJ_yr => 0.0,
              BB_radiation_at_Temp_in_atm_ZJ_yr => 0.0,
              Blocked_by_CH4 => 0.0,
              Blocked_by_CO2 => 0.0,
              Blocked_by_H20 => 0.0,
              Blocked_by_H20_future_linear_equ => 0.0,
              Blocked_by_H20_future_poly_equ => 0.0,
              Blocked_by_H20_future_poly_equ_dyn => 0.0,
              Blocked_by_H20_future_poly_equ_dyn_0 => 0.0,
              Blocked_by_H20_hist_Table_lookup => 0.0,
              Blocked_by_h20_poly_used => 0.0,
              Blocked_by_H20_Table_lookup => 0.0,
              Blocked_by_H2O_hist_and_fut => 0.0,
              Blocked_by_H2O_poly_dyn => 0.0,
              Blocked_by_H2O_poly_equ => 0.0,
              Blocked_by_otherGHG => 0.0,
              Blocking_multiplier_from_Kyoto_Flour => 0.0,
              Blocking_multiplier_from_Montreal_gases => 0.0,
              Blocking_multiplier_from_N2O => 0.0,
              Blocking_of_LW_rad_by_clouds => 0.0,
              C_absorption_by_ocean_from_atm_for_accumulation => 0.0,
              C_diffusion_into_ocean_from_atm => 0.0,
              C_diffusion_into_ocean_from_atm_MtC_yr => 0.0,
              C_in_biomass => 0.0,
              C_in_GRASS_DeadB_and_soil_GtC => 0.0,
              C_in_GRASS_GtC => 0.0,
              C_in_GRASS_LB_GtC => 0.0,
              C_in_NF_DeadB_and_soil_GtC => 0.0,
              C_in_NF_GtC => 0.0,
              C_in_NF_LB_GtC => 0.0,
              C_in_TROP_DeadB_and_soil_GtC => 0.0,
              C_in_TROP_GtC => 0.0,
              C_in_TROP_LB_GtC => 0.0,
            ]
          end
          push!(startEquationConstructors, generateStartEquations9)
        end
        begin
          function generateStartEquations10()
            [
              C_in_TUNDRA_DeadB_and_soil_GtC => 0.0,
              C_in_TUNDRA_GtC => 0.0,
              C_in_TUNDRA_LB_GtC => 0.0,
              C_release_from_permafrost_melting_as_CO2_GtC_yr => 0.0,
              C_released_from_permafrost_as_either_CH4_or_CO2_in_GtC_yr => 0.0,
              C_removal_rate_from_atm_for_nature_May_2020_GtC_y => 0.0,
              C_runoff_from_biomass_soil => 0.0,
              Carbon_captured_and_stored_GtC___yr => 0.0,
              Carbon_concentration_in_cold_surface_ocean => 0.0,
              Carbon_concentration_in_CWTtB => 0.0,
              Carbon_concentration_in_deep_box_GtC_per_G_cubicM => 0.0,
              Carbon_concentration_in_intermdiate_box_GtC_per_G_cubicM => 0.0,
              Carbon_concentration_in_warm_surface => 0.0,
              Carbon_flow_from_cold_surface_downwelling_Gtc_per_yr => 0.0,
              Carbon_flow_from_cold_to_deep_GtC_per_yr => 0.0,
              Carbon_flow_from_deep => 0.0,
              Carbon_flow_from_intermediate_to_surface_box_GtC_per_yr => 0.0,
              Carbon_flow_from_warm_to_cold_surface_GtC_per_yr => 0.0,
              Carbon_in_cold_ocean_0_to_100m_1850_GtC => 0.0,
              Carbon_in_cold_ocean_trunk_100m_to_bottom_1850_GtC => 0.0,
              Carbon_in_ocean_deep_1k_to_bottom_ocean_1850_GtC => 0.0,
              Carbon_in_ocean_upwelling_100m_to_1km_1850_GtC => 0.0,
              Carbon_in_top_ocean_layer_1850_GtC => 0.0,
              Carbon_in_top_ocean_layer_GtC => 0.0,
              Carbon_in_warm_ocean_0_to_100m_1850_GtC => 0.0,
              CC_in_cold_downwelling_ymoles_per_litre => 0.0,
              CC_in_cold_downwelling_ymoles_per_litre__dimensionless_ => 0.0,
              CC_in_cold_surface_ymoles_per_litre => 0.0,
              CC_in_cold_surface_ymoles_per_litre__dimensionless_ => 0.0,
              CC_in_deep_box_ymoles_per_litre => 0.0,
              CC_in_deep_box_ymoles_per_litre__dimensionless_ => 0.0,
              CC_in_intermediate_box_ymoles_per_litre => 0.0,
              CC_in_intermediate_box_ymoles_per_litre__dimensionless_ => 0.0,
              CC_in_warm_surface_ymoles_per_litre => 0.0,
              CC_in_warm_surface_ymoles_per_litre__dimensionless_ => 0.0,
              CH4_all_emissions_GtC_yr => 0.0,
              CH4_concentration_ppb => 0.0,
              CH4_conversion_to_CO2_and_H2O => 0.0,
              CH4_emissions_before_co2e_exp => 0.0,
              CH4_emissions_CO2e_after_exp => 0.0,
              CH4_emissions_CO2e_after_exp_12a => 0.0,
              CH4_emissions_from_wetlands_destruction => 0.0,
              CH4_in_permafrost_area_melted___frozen_before_heat_constraint_GtC_yr => 0.0,
              CH4_in_the_atmosphere_converted_to_CO2 => 0.0,
              CH4_per_sqkm_of_wetlands => 0.0,
              CH4_release___capture_from_permafrost_area_loss___gain_GtC_yr => 0.0,
              CO2_conc_atm_less_CO2_conc_sea => 0.0,
              CO2_conc_in_cold_surface_water_in_ppm => 0.0,
              CO2_conc_in_warm_surface_water_in_ppm => 0.0,
              CO2_concentration_calculated_as_a_1pct_pa_exponential_increase_ppm => 0.0,
            ]
          end
          push!(startEquationConstructors, generateStartEquations10)
        end
        begin
          function generateStartEquations11()
            [
              CO2_concentration_ppm => 0.0,
              CO2_concentration_used__after_any_experiments__ppm => 0.0,
              CO2_emissions_before_co2e_exp => 0.0,
              CO2_emissions_CO2e_after_exp => 0.0,
              CO2_flow_from_GRASS_to_atmosphere_GtC_yr => 0.0,
              CO2_flow_from_NF_to_atmosphere_GtC_yr => 0.0,
              CO2_flow_from_TROP_to_atmosphere_GtC_yr => 0.0,
              CO2_flow_from_TUNDRA_to_atmosphere_GtC_yr => 0.0,
              CO2_flux_from_atm_to_GRASS_for_new_growth_GtC_yr => 0.0,
              CO2_flux_from_atm_to_NF_for_new_growth_GtC_yr => 0.0,
              CO2_flux_from_atm_to_TROP_for_new_growth_GtC_yr => 0.0,
              CO2_flux_from_atm_to_TUNDRA_for_new_growth_GtC_yr => 0.0,
              CO2_flux_GRASS_to_atm_Gtc_yr => 0.0,
              CO2_flux_NF_to_atm_Gtc_yr => 0.0,
              CO2_flux_TROP_to_atm_GtC_yr => 0.0,
              CO2_flux_TUNDRA_to_atm_Gtc_yr => 0.0,
              CO2_radiative_forcing_since_1850_using_Myhre_formula_W_pr_m2 => 0.0,
              Cold_dense_water_sinking_in_Sverdrup => 0.0,
              Concentration_of_C_in_ocean_top_layer_in_1850 => 0.0,
              Contrib_of_BARREN_land_to_albedo_land => 0.0,
              Contrib_of_GRASS_to_albedo_land => 0.0,
              Contrib_of_ICE_ON_LAND_to_albedo_land => 0.0,
              Contrib_of_NF_to_albedo_land => 0.0,
              Contrib_of_TROP_to_albedo_land => 0.0,
              Contrib_of_TUNDRA_to_albedo_land => 0.0,
              Contribution_to_forcing_by_CH4 => 0.0,
              Contribution_to_forcing_by_CO2 => 0.0,
              Contribution_to_forcing_by_H2O => 0.0,
              Contribution_to_forcing_by_othGHG => 0.0,
              Convection_aka_sensible_heat_flow => 0.0,
              Convection_aka_sensible_heat_flow_W_m2 => 0.0,
              Convection_as_f_of_temp_ZJ_yr => 0.0,
              Conversion_constant_GtC_to_ppm => 0.0,
              Conversion_constant_heat_ocean_deep_to_temp => 0.0,
              Conversion_heat_atm_to_temp => 0.0,
              Conversion_heat_surface_to_temp => 0.0,
              dbl_CO2_exp => 0.0,
              Deep_ocean__cold__volume => 0.0,
              delta_C_in_atmosphere_GtC => 0.0,
              delta_C_in_biomass_GtC => 0.0,
              delta_C_in_ocean_GtC => 0.0,
              delta_CO2_concentration_since_1850_ppm => 0.0,
              delta_Temp_deep_ocean_degC => 0.0,
              Depositing_of_C_to_sediment => 0.0,
              Effect_of_acidification_on_NMPP => 0.0,
              Effect_of_C_concentration_on_NMPP => 0.0,
              Effect_of_CO2_on_new_biomass_growth => 0.0,
              Effect_of_heat_in_atm_on_melting_ice__cut_off_ => 0.0,
              Effect_of_humidity_on_shifting_biomes => 0.0,
              Effect_of_population_and_urbanization_on_biomass_use => 0.0,
            ]
          end
          push!(startEquationConstructors, generateStartEquations11)
        end
        begin
          function generateStartEquations12()
            [
              Effect_of_temp_on_melting_antarctic_ice => 0.0,
              Effect_of_temp_on_melting_greenland_ice => 0.0,
              Effect_of_temp_on_melting_greenland_ice_that_slid_into_the_ocean => 0.0,
              Effect_of_temp_on_melting_or_freezing_glacial_ice => 0.0,
              Effect_of_temp_on_melting_or_freezing_of_Arctic_ice => 0.0,
              Effect_of_temperature_on_fire_incidence_dimensionless => 0.0,
              Effect_of_temperature_on_new_biomass_growth_dimensionless => 0.0,
              Effect_of_temperature_on_NMPP => 0.0,
              Effective_time_to_melt_Arctic_ice_at_the_reference_delta_temp => 0.0,
              Effective_time_to_melt_glacial_ice_at_the_reference_delta_temp => 0.0,
              Effective_time_to_melt_greenland_ice_at_the_reference_delta_temp => 0.0,
              Effective_time_to_melt_or_freeze_antarctic_ice_at_the_reference_delta_temp => 0.0,
              Effective_Time_to_regrow_TROP_after_deforesting_yr => 0.0,
              Emissions_of_aerosols_1850_to_2100_with_IPCC_Fig_pg_1037_Exp => 0.0,
              Emissions_of_anthro_CH4_1850_to_2100_GtC_yr => 0.0,
              Emissions_of_anthro_CH4_1850_to_2100_linearly_reduced_from_2015_GtC_yr => 0.0,
              Emissions_of_anthro_CH4_1850_to_2100_RCP_with_JR_shape_forecast_GtC_yr => 0.0,
              Emissions_of_anthro_CH4_1850_to_2100_RCP_with_Xpct_annual_incr_forecast_GtC_yr => 0.0,
              Emissions_of_anthro_CH4_1850_to_2100_RCP3_GtC_yr => 0.0,
              Emissions_of_anthro_CH4_1850_to_2100_RCP45_GtC_yr => 0.0,
              Emissions_of_anthro_CH4_1850_to_2100_RCP6_GtC_yr => 0.0,
              Emissions_of_anthro_CH4_1850_to_2100_RCP85_GtC_yr => 0.0,
              Emissions_of_anthro_CO2_1850_to_2100_linearly_reduced_from_2015_GtC_yr => 0.0,
              Emissions_of_anthro_CO2_1850_to_2100_RCP_with_Xpct_annual_incr_forecast_GtC_yr => 0.0,
              Emissions_of_CH4_1850_to_2100_GtC_yr_with_IPCC_Fig_pg_1037_Exp => 0.0,
              Emissions_of_CO2_1850_to_2100_GtC_yr_with_EXP_12a => 0.0,
              Emissions_of_CO2_1850_to_2100_GtC_yr_with_IPCC_Fig_pg_1037_Exp => 0.0,
              Emissions_of_CO2_1850_to_2100_GtC_yr => 0.0,
              Evaporation_aka_latent_heat_flow => 0.0,
              Evaporation_aka_latent_heat_flow_W_m2 => 0.0,
              Evaporation_as_f_of_temp_ZJ_yr => 0.0,
              Exogenous_sliding_of_Greenland_ice_into_the_ocean => 0.0,
              Exp_12a_reduction_in_emissions => 0.0,
              Exp_12a_reduction_in_emissions_LOOKUP => 0.0,
              EXP_12b_CCS_from_2015 => 0.0,
              EXP_12c_stopping_TROP_deforestation_from_2015 => 0.0,
              EXP_12e_white_surfaces_ease_in => 0.0,
              exp0 => 0.0,
              exp0_dyn => 0.0,
              exp1 => 0.0,
              exp1_dyn => 0.0,
              exp2 => 0.0,
              exp2_dyn => 0.0,
              exp3 => 0.0,
              exp3_dyn => 0.0,
              Experimental_doubling_of_constant_C_emissions => 0.0,
              Experimental_release_of_constant_fossil_C_emissions_GtC_yr => 0.0,
              Experimental_release_of_methane => 0.0,
              f_M_1750_N_2010__for_ch4_forcing => 0.0,
              f_M_2010_N_cur_ => 0.0,
            ]
          end
          push!(startEquationConstructors, generateStartEquations12)
        end
        begin
          function generateStartEquations13()
            [
              f_M_cur_N_2010_ => 0.0,
              f_M2010_N_1750__for_n20_forcing => 0.0,
              Flow_from_atm_to_biomass_GtC_pr_yr => 0.0,
              Flow_from_biomass_to_atm_Gtc_pr_yr => 0.0,
              Flow_of_cold_surface_water_welling_down_GcubicM_per_yr => 0.0,
              Flow_of_cold_water_sinking_to_very_bottom_GcubicM_per_yr => 0.0,
              Flow_of_heat_to_atm_ZJ_yr => 0.0,
              Flow_of_heat_to_deep_ocean => 0.0,
              Flow_of_heat_to_deep_ocean_btw_72_and_08 => 0.0,
              Flow_of_heat_to_surface_ocean => 0.0,
              Flow_of_heat_to_surface_ocean_btw_1972_and_2008 => 0.0,
              Flow_of_water_from_warm_to_cold_surface_Gcubicm_per_yr => 0.0,
              for_display_yr_on_yr_change_in_C_in_ocean_GtC_yr => 0.0,
              Frac_atm_absorption => 0.0,
              Frac_blocked_by_ALL_GHG => 0.0,
              Frac_blocked_by_ALL_GHG_LESS_watervapor => 0.0,
              Frac_vol_cold_ocean_0_to_100m_of_total => 0.0,
              Frac_vol_cold_ocean_downwelling_of_total => 0.0,
              Frac_vol_deep_ocean_of_total => 0.0,
              Frac_vol_ocean_upwelling_of_total => 0.0,
              Frac_vol_warm_ocean_0_to_100m_of_total => 0.0,
              Fraction_blocked_by_CH4_spectrum => 0.0,
              Fraction_blocked_by_CO2_spectrum => 0.0,
              Fraction_blocked_by_other_GHG => 0.0,
              Fraction_GRASS_being_deforested_1_yr => 0.0,
              Fraction_of_C_released_from_permafrost_released_as_CH4_dmnl => 0.0,
              Fraction_of_ocean_classified_as_cold_surface => 0.0,
              Fraction_TUNDRA_being_deforested_1_yr => 0.0,
              Future_shape_of_anthropogenic_aerosol_emissions => 0.0,
              Ga__BB_radiation_less_TOA_radiation_W_m2 => 0.0,
              Glacial_ice_area_decrease_Mkm2_pr_yr => 0.0,
              Glacial_ice_area_increase_Mkm2_pr_yr => 0.0,
              Glacial_ice_area_km2 => 0.0,
              Glacial_ice_freezing_km3_yr => 0.0,
              Glacial_ice_melting__pos__or_freezing__neg__km3_yr => 0.0,
              Glacial_ice_melting_as_water_km3_yr => 0.0,
              Glacial_ice_melting_km3_yr => 0.0,
              GRASS_being_deforested_Mkm2_yr => 0.0,
              GRASS_being_harvested_Mkm2_yr => 0.0,
              GRASS_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr => 0.0,
              GRASS_Biomass_in_construction_material_being_burnt_GtBiomass_yr => 0.0,
              GRASS_Biomass_in_construction_material_left_to_rot_GtBiomass_yr => 0.0,
              GRASS_biomass_new_growing_GtBiomass___yr => 0.0,
              GRASS_burning_Mkm2_yr => 0.0,
              GRASS_Dead_biomass_decomposing_GtBiomass_yr => 0.0,
              GRASS_DeadB_and_SOM_tB_per_km2 => 0.0,
              GRASS_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr => 0.0,
              GRASS_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr => 0.0,
              GRASS_for_construction_use_GtBiomass_yr => 0.0,
              GRASS_historical_deforestation_pct_yr => 0.0,
            ]
          end
          push!(startEquationConstructors, generateStartEquations13)
        end
        begin
          function generateStartEquations14()
            [
              GRASS_land_taken_out_of_use_GtBiomass => 0.0,
              GRASS_land_taken_out_of_use_Mkm2 => 0.0,
              GRASS_living_biomass_densitiy_tBiomass_pr_km2 => 0.0,
              GRASS_Living_biomass_rotting_GtBiomass_yr => 0.0,
              GRASS_potential_less_actual_living_biomass_GtBiomass => 0.0,
              GRASS_potential_living_biomass_GtBiomass => 0.0,
              GRASS_regrowing_after_being_burnt_Mkm2_yr => 0.0,
              GRASS_regrowing_after_being_deforested_Mkm2_yr => 0.0,
              GRASS_regrowing_after_harvesting_Mkm2_yr => 0.0,
              GRASS_runoff => 0.0,
              GRASS_soil_degradation_from_forest_fires_GtBiomass_yr => 0.0,
              GRASS_with_normal_cover_Mkm2 => 0.0,
              Greenland_ice_area_decrease_Mkm2_pr_yr => 0.0,
              Greenland_ice_area_increase_Mkm2_pr_yr => 0.0,
              Greenland_ice_area_km2 => 0.0,
              Greenland_ice_freezing_km3_yr => 0.0,
              Greenland_ice_losing__pos__or_gaining__neg__GtIce_yr => 0.0,
              Greenland_ice_melting__pos__or_freezing__neg__km3_yr => 0.0,
              Greenland_ice_melting_as_water_km3_yr => 0.0,
              Greenland_ice_melting_km3_yr => 0.0,
              Greenland_ice_melting_that_slid_into_the_ocean_km3_yr => 0.0,
              Greenland_ice_sliding_into_the_ocean_km3_yr => 0.0,
              Greenland_ice_that_slid_into_the_ocean_melting__pos__or_freezing__neg__km3_yr => 0.0,
              Guldberg_Waage_air_sea_formulation => 0.0,
              Heat_actually_gained___needed_for_freezing___unfreezing_of_permafrost_ZJ_yr => 0.0,
              Heat_flow_from_the_earths_core => 0.0,
              Heat_gained___needed_for_the_desired_freezing___unfreezing_of_permafrost_ZJ_yr => 0.0,
              Heat_used_in_melting__pos__or_freezing__neg__antarctic_ice_ZJ_yr => 0.0,
              Heat_used_in_melting__pos__or_freezing__neg__arctic_sea_ice_ZJ_yr => 0.0,
              Heat_used_in_melting__pos__or_freezing__neg__glacial_ice_ZJ_yr => 0.0,
              Heat_used_in_melting__pos__or_freezing__neg__Greenland_ice_that_slid_into_the_water_ZJ_yr => 0.0,
              Heat_used_in_melting__pos__or_freezing__neg__Greenland_ice_ZJ_yr => 0.0,
              Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_antarctic_ice_W_m2 => 0.0,
              Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_antarctic_ice_ZJ_yr => 0.0,
              Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_arctic_ice_W_m2 => 0.0,
              Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_arctic_ice_ZJ_yr => 0.0,
              Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_Greenland_ice_that_slid_into_the_ocean_W_m2 => 0.0,
              Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_Greenland_ice_that_slid_into_the_ocean_ZJ_yr => 0.0,
              Heat_withdrawn_from_atm_by_melting__pos____added__neg__by_freezing_ice_W_m2 => 0.0,
              Heat_withdrawn_from_atm_by_melting__pos____added__neg__by_freezing_ice_ZJ_yr => 0.0,
              HI_clouds_net_effect__pos_warming__neg_cooling__W_m2 => 0.0,
              Hist_Frac_atm_absorption => 0.0,
              Human_activity_CH4_emissions => 0.0,
              Human_activity_CH4_emissions_GtCO2e_yr => 0.0,
              Humidity_of_atmosphere_current_g_kg => 0.0,
              Humidity_of_atmosphere_g_kg => 0.0,
              Ice_on_land_area_Mkm2 => 0.0,
              Incoming_solar_W_m2 => 0.0,
              Incoming_solar_ZJ_yr => 0.0,
              InputEmissions_for_tipping_point_search => 0.0,
            ]
          end
          push!(startEquationConstructors, generateStartEquations14)
        end
        begin
          function generateStartEquations15()
            [
              Intercept_blocked_by_H20_future_equ => 0.0,
              Kyoto_Flour_concentration_ppt => 0.0,
              Kyoto_Flour_degradation => 0.0,
              Kyoto_Flour_emissions => 0.0,
              Kyoto_Flour_emissions_after_exp => 0.0,
              Kyoto_Flour_emissions_after_exp_12a => 0.0,
              Kyoto_Flour_emissions_before_exp => 0.0,
              Kyoto_Flour_emissions_GtCO2e_yr => 0.0,
              Kyoto_Flour_emissions_RCPs_or_JR52 => 0.0,
              Land_area_km2 => 0.0,
              Land_covered_with_ice_km2 => 0.0,
              Land_covered_with_ice_Mkm2 => 0.0,
              LO_clouds_net_effect__pos_warming__neg_cooling__W_m2 => 0.0,
              LW_Blocking_multiplier_from_other_GHG => 0.0,
              LW_Clear_sky_emissions_from_atm => 0.0,
              LW_Clear_sky_emissions_from_atm_W_m2 => 0.0,
              LW_clear_sky_emissions_to_surface => 0.0,
              LW_clear_sky_emissions_to_surface_W_m2 => 0.0,
              LW_Cloudy_sky_emissions_from_atm => 0.0,
              LW_Cloudy_sky_emissions_from_atm_W_m2 => 0.0,
              LW_HI_cloud_radiation => 0.0,
              LW_HI_cloud_radiation_reference_in_1850_W_m2 => 0.0,
              LW_HI_cloud_radiation_W_m2 => 0.0,
              LW_LO_cloud_radiation => 0.0,
              LW_LO_cloud_radiation_W_m2 => 0.0,
              LW_radiation_blocked_by_CH4__pct_ => 0.0,
              LW_radiation_blocked_by_CO2__pct_ => 0.0,
              LW_radiation_blocked_by_H2O__pct_ => 0.0,
              LW_radiation_blocked_by_other_GHG__pct_ => 0.0,
              LW_re_radiated_by_clouds => 0.0,
              LW_re_radiated_by_clouds_W_m2 => 0.0,
              LW_surface_emission => 0.0,
              LW_surface_emission_W_m2 => 0.0,
              LW_surface_emissions_escaping_through_atm_window => 0.0,
              LW_surface_emissions_NOT_escaping_through_atm_window => 0.0,
              LW_surface_emissions_NOT_escaping_through_atm_window_W_m2 => 0.0,
              LW_TOA_radiation_from_atm_to_space => 0.0,
              LW_TOA_radiation_from_atm_to_space_difference_wrt_1850 => 0.0,
              LW_TOA_radiation_from_atm_to_space_difference_wrt_1850_W_m2 => 0.0,
              LW_TOA_radiation_from_atm_to_space_W_m2 => 0.0,
              M_2010 => 0.0,
              M_cur => 0.0,
              Man_made_CH4_emissions_pct => 0.0,
              Man_made_fossil_C_emissions_for_cumulation_GtC_yr => 0.0,
              Man_made_fossil_C_emissions_GtC_yr => 0.0,
              Man_made_fossil_C_emissions_GtCO2e_yr => 0.0,
              Melting_constraint_from_the_heat_in__ocean__surface_reservoir => 0.0,
              Melting_constraint_from_the_heat_in_atmosphere_reservoir_fraction => 0.0,
              Melting_restraint_for_permafrost_from_heat_in_atmophere => 0.0,
              Methane_hydrates_released_and_converted_to_CO2_by_bacteria_GtC => 0.0,
            ]
          end
          push!(startEquationConstructors, generateStartEquations15)
        end
        begin
          function generateStartEquations16()
            [
              Methanehydrate_experimental_release_GtC__yr => 0.0,
              MODEL_CH4_in_atm_in_ppb => 0.0,
              MODEL_CO2_concentration_in_atmosphere2_ppm => 0.0,
              Model_Volcanic_aerosol_forcing_W_m2 => 0.0,
              Montreal_emissions_GtCO2e_yr => 0.0,
              Montreal_gases_concentration_ppt => 0.0,
              Montreal_gases_degradation => 0.0,
              Montreal_gases_emissions => 0.0,
              Montreal_gases_emissions_after_exp_12a => 0.0,
              Montreal_gases_emissions_before_exp => 0.0,
              Montreal_gases_emissions_CO2e_after_exp => 0.0,
              Montreal_gases_emissions_RCPs_or_JR52 => 0.0,
              N_2010 => 0.0,
              N_cur => 0.0,
              N20_emissions_RCPs_or_JR52 => 0.0,
              N2O_concentration_ppb => 0.0,
              N2O_degradation_MtN2O_yr => 0.0,
              N2O_man_made_emissions => 0.0,
              N2O_man_made_emissions_after_exp => 0.0,
              N2O_man_made_emissions_exp_12a => 0.0,
              N2O_man_made_emissions_GtCO2e_yr => 0.0,
              NatEvent_d__slowing_down_ocean_circulation_from_2015 => 0.0,
              Natural_CH4_emissions => 0.0,
              Natural_CH4_emissions_pct => 0.0,
              NATURE_CCS_Fig3_GtC_yr => 0.0,
              NATURE_CCS_removal_experiment_multiplier => 0.0,
              Net_additions_to_C_in_TUNDRA_DeadB_and_soil_GtC => 0.0,
              Net_additions_to_C_in_TUNDRA_LB_GtC => 0.0,
              Net_C_flow_from_atm_to_biomass_GtC_pr_yr => 0.0,
              Net_C_to_atm => 0.0,
              Net_C_to_atm_rate => 0.0,
              Net_CO2_flow_between_grass_and_atmosphere_GtC => 0.0,
              Net_CO2_flow_between_TUNDRA_and_atmosphere_GtC => 0.0,
              Net_flow_of_heat_into_surface => 0.0,
              Net_flux_to_ocean_GtC_yr => 0.0,
              Net_heat_flow_ocean_between_surface_and_deep_per_K_of_difference_ZJ_yr_K => 0.0,
              Net_heat_flow_ocean_from_surface_to_deep__ZJ_yr_ => 0.0,
              Net_heat_flow_ocean_from_surface_to_deep_W_m2 => 0.0,
              Net_heat_flow_to_atm_ZJ_yr__needed_for_comparisons_with_history_ => 0.0,
              Net_marine_primary_production_NMPP_GtC_pr_yr => 0.0,
              NEW_Temp_ocean_surface_in_1850_in_K => 0.0,
              NF_Avg_life_biomass_yr => 0.0,
              NF_being_deforested_Mkm2_yr => 0.0,
              NF_being_harvested_by_clear_cutting_Mkm2_yr => 0.0,
              NF_being_harvested_Mkm2_yr => 0.0,
              NF_being_harvested_normally_Mkm2_yr => 0.0,
              NF_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr => 0.0,
              NF_Biomass_in_construction_material_being_burnt_GtBiomass_yr => 0.0,
              NF_biomass_new_growing_GtBiomass___yr => 0.0,
              NF_burning_Mkm2_yr => 0.0,
            ]
          end
          push!(startEquationConstructors, generateStartEquations16)
        end
        begin
          function generateStartEquations17()
            [
              NF_clear_cut_fraction => 0.0,
              NF_Dead_biomass_decomposing_GtBiomass_yr => 0.0,
              NF_DeadB_and_SOM_tB_per_km2 => 0.0,
              NF_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr => 0.0,
              NF_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr => 0.0,
              NF_for_construction_use_GtBiomass_yr => 0.0,
              NF_historical_deforestation_pct_yr => 0.0,
              NF_land_taken_out_of_use_GtBiomass => 0.0,
              NF_land_taken_out_of_use_Mkm2 => 0.0,
              NF_living_biomass_densitiy_tBiomass_pr_km2 => 0.0,
              NF_Living_biomass_rotting_GtBiomass_yr => 0.0,
              NF_potential_less_actual_living_biomass_GtBiomass => 0.0,
              NF_potential_living_biomass_GtBiomass => 0.0,
              NF_regrowing_after_being_burnt_Mkm2_yr => 0.0,
              NF_regrowing_after_being_clear_cut_Mkm2_yr => 0.0,
              NF_regrowing_after_being_deforested_Mkm2_yr => 0.0,
              NF_regrowing_after_harvesting_Mkm2_yr => 0.0,
              NF_runoff => 0.0,
              NF_soil_degradation_from_clear_cutting_GtBiomass_yr => 0.0,
              NF_soil_degradation_from_forest_fires_GtBiomass_yr => 0.0,
              NF_Speed_of_regrowth_yr => 0.0,
              NF_Sum_outflows_GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass_yr => 0.0,
              NF_TUNDRA_Biomass_in_construction_material_left_to_rot_GtBiomass_yr => 0.0,
              NF_usage_as_pct_of_potial_area => 0.0,
              NF_usage_cutoff => 0.0,
              NF_with_normal_cover_Mkm2 => 0.0,
              Ocean_area_km2 => 0.0,
              Ocean_circulation_slowdown_from_Greenland_ice_sliding_into_the_Atlantic => 0.0,
              Ocean_heat_used_for_melting_ZJ_yr => 0.0,
              Ocean_surface_area_km2 => 0.0,
              Ocean_surface_delta_temp_to_1850_C => 0.0,
              Open_water_as_frac_of_ocean_area => 0.0,
              Outgoing_radiation_at_TOA_W_m2 => 0.0,
              pct_change_in_fraction_blocked_by_ALL_GHG_wrt_1850 => 0.0,
              pct_change_in_fraction_blocked_by_C02_wrt_1850 => 0.0,
              pct_change_in_fraction_blocked_by_CH4_wrt_1850 => 0.0,
              pct_change_in_fraction_blocked_by_othGHG_wrt_1850 => 0.0,
              pct_reduction_in_C_in_GRASS => 0.0,
              pct_reduction_in_C_in_NF => 0.0,
              pct_reduction_in_C_in_TROP => 0.0,
              pct_reduction_in_C_in_TUNDRA => 0.0,
              Permafrost_area_km2 => 0.0,
              Permafrost_CH4_emissions_pct => 0.0,
              Permafrost_melting_cutoff => 0.0,
              pH_in_cold_deep_water => 0.0,
              ph_in_cold_downwelling_water => 0.0,
              pH_in_cold_suface_water => 0.0,
              pH_in_surface => 0.0,
              pH_in_upwelling_water => 0.0,
              pH_in_warm_surface_water => 0.0,
            ]
          end
          push!(startEquationConstructors, generateStartEquations17)
        end
        begin
          function generateStartEquations18()
            [
              POLICY_4_Stopping_logging_in_Northern_forests => 0.0,
              Radiation_balance_at_TOA_in_less_out_W_m2 => 0.0,
              Radiative_forcing_from_CH4_wrt_1850_W_m2 => 0.0,
              Radiative_forcing_from_CO2_wrt_1850_W_m2 => 0.0,
              Radiative_forcing_from_H2O_wrt_1850_W_m2 => 0.0,
              Radiative_forcing_from_othGHG_wrt_1850_W_m2 => 0.0,
              Radiative_forcing_wrt_1850_W_m2_0 => 0.0,
              Rate_of_destruction_of_wetlands => 0.0,
              Ratio_of_area_covered_by_high_clouds_current_to_1850 => 0.0,
              Ratio_of_area_covered_by_low_clouds_current_to_1850 => 0.0,
              RCPFossil_fuel_usage_cutoff => 0.0,
              Reflected_Solar_SW => 0.0,
              Reflected_Solar_SW_W_m2 => 0.0,
              RF_CH4_IPCC_formula_W_m2 => 0.0,
              RF_CO2_Model_Myhre_formula => 0.0,
              RF_CO2_Model_Myhre_formula_1850 => 0.0,
              RF_CO2_RCP3_Myhre_formula => 0.0,
              RF_CO2_RCP45_Myhre_formula => 0.0,
              RF_CO2_RCP6_Myhre_formula => 0.0,
              RF_CO2_RCP85_Myhre_formula => 0.0,
              RF_from_Ga__BB_radiation_less_TOA_radiation_W_m2 => 0.0,
              RF_N20_IPCC_formula_W_m2 => 0.0,
              Sea_level_change_from_melting_ice_and_thermal_expansion_m => 0.0,
              Sea_level_change_from_thermal_expansion_deep_m => 0.0,
              Sea_level_change_from_thermal_expansion_surface_m => 0.0,
              Sea_level_rise_from_melting_ice_m => 0.0,
              Sea_level_rise_history_m => 0.0,
              Seconds_per_yr => 0.0,
              Sensitivity_of_high_cloud_coverage_to_temp => 0.0,
              Shifting_GRASS_to_DESERT_Mkm2_yr => 0.0,
              Shifting_GRASS_to_NF_Mkm2_yr => 0.0,
              Shifting_GRASS_to_TROP_Mkm2_yr => 0.0,
              Shifting_ice_on_land_to_tundra_Mkm2_yr => 0.0,
              Shifting_ice_to_tundra_from_detail_ice_on_land_Mkm2_pr_yr => 0.0,
              Shifting_NF_to_GRASS_Mkm2_yr => 0.0,
              Shifting_NF_to_TROP_Mkm2_yr => 0.0,
              Shifting_NF_to_Tundra_Mkm2_yr => 0.0,
              Shifting_TROP_to_GRASS_Mkm2_yr => 0.0,
              Shifting_TROP_to_NF_Mkm2_yr => 0.0,
              Shifting_tundra_to_ice_from_detail_ice_on_land_Mkm2_pr_yr => 0.0,
              Shifting_tundra_to_ice_on_land_Mkm2_yr => 0.0,
              Shifting_Tundra_to_NF_Mkm2_yr => 0.0,
              SHUT_OFF_permafrost => 0.0,
              Sifting_DESERT_to_GRASS_Mkm2_yr => 0.0,
              Slider_for_H2O_slope => 0.0,
              Slope_blocked_by_H20_future_equ => 0.0,
              Slope_btw_temp_and_permafrost_melting___freezing => 0.0,
              Slope_of_effect_of_temp_shifting_DESERT_to_GRASS => 0.0,
              Slope_temp_vs_glacial_ice_melting => 0.0,
              Slowing_of_recapture_of_CH4_dmnl => 0.0,
            ]
          end
          push!(startEquationConstructors, generateStartEquations18)
        end
        begin
          function generateStartEquations19()
            [
              Snowball_earth_cutoff => 0.0,
              Solar_cycle_W_m2 => 0.0,
              Solar_sine_forcing_W_m2 => 0.0,
              Stop_of_human_deforestation => 0.0,
              Sum_biomes_Mkm2 => 0.0,
              sum_blocked => 0.0,
              Sum_heat_to_ocean_1972_to_2008_ZJ => 0.0,
              Sum_outflows_GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass_yr => 0.0,
              Surface_deep__ocean__temp_diff_degC => 0.0,
              Surface_imbalance_pos_is_TO_surface => 0.0,
              Surface_imbalance_pos_is_TO_surface_W_m2 => 0.0,
              Surface_ocean__warm__volume => 0.0,
              SW_Atmospheric_absorption => 0.0,
              SW_Atmospheric_absorption_W_m2 => 0.0,
              SW_clear_sky_reflection_aka_scattering => 0.0,
              SW_clear_sky_reflection_aka_scattering_W_m2 => 0.0,
              SW_HI_cloud_efffect_aka_cloud_albedo => 0.0,
              SW_HI_cloud_efffect_aka_TOA_albedo_W_m2 => 0.0,
              SW_LO_cloud_efffect_aka_cloud_albedo => 0.0,
              SW_LO_cloud_efffect_aka_cloud_albedo_W_m2 => 0.0,
              SW_surface_absorption => 0.0,
              SW_surface_absorption_W_m2_wrt_1850 => 0.0,
              SW_surface_absorption_W_m2 => 0.0,
              SW_surface_reflection => 0.0,
              SW_surface_reflection_W_m2_wrt_1850 => 0.0,
              SW_surface_reflection_W_m2 => 0.0,
              SW_to_surface => 0.0,
              SW_to_surface_W_m2 => 0.0,
              Temp__ocean__deep_in_1850_in_K => 0.0,
              Temp__ocean__deep_in_C => 0.0,
              Temp__ocean__surface_in_K => 0.0,
              Temp_atm_average_K => 0.0,
              Temp_atm_in_C => 0.0,
              Temp_driver_to_shift_biomes_degC => 0.0,
              Temp_gradient => 0.0,
              Temp_gradient_minus_1 => 0.0,
              Temp_gradient_minus_1___slope => 0.0,
              Temp_ocean_deep_in_K => 0.0,
              Temp_of_cold_downwelling_water => 0.0,
              Temp_of_cold_surface_water => 0.0,
              Temp_surface_anomaly_compared_to_1850_degC => 0.0,
              Temp_surface_average_K => 0.0,
              Temp_surface_C => 0.0,
              Temp_surface_current_divided_by_value_in_1850_K_K => 0.0,
              Thermal_expansion_deep_in_1850_pct => 0.0,
              Thermal_expansion_deep_pct => 0.0,
              Thermal_expansion_surface_in_1850_pct => 0.0,
              Thermal_expansion_surface_pct => 0.0,
              Time_in_trunk => 0.0,
              Time_less_Greenland_slide_experiment_start_yr => 0.0,
            ]
          end
          push!(startEquationConstructors, generateStartEquations19)
        end
        begin
          function generateStartEquations20()
            [
              Time_to_degrade_Kyoto_Flour_yr => 0.0,
              Time_to_regrow_NF_after_buning_yr => 0.0,
              Tipping_point_search_emissions_GtCO2e_yr => 0.0,
              Tipping_point_year_of_peak => 0.0,
              Total_carbon_in_Ocean_1850_GtC => 0.0,
              Total_carbon_in_ocean_GtC => 0.0,
              Total_CO2e_emissions_as_f_peak__GtCO2e_yr => 0.0,
              Total_net_aerosol_forcing_ZJ_yr => 0.0,
              Total_net_aerosol_forcings_W_m2 => 0.0,
              Total_sea_level_change_from_thermal_expansion_m => 0.0,
              Total_volume_of_ocean_water_GcubicM => 0.0,
              TROP_being_deforested_Mkm2_yr => 0.0,
              TROP_being_harvested_by_clear_cutting_Mkm2_yr => 0.0,
              TROP_being_harvested_Mkm2_yr => 0.0,
              TROP_being_harvested_normally_Mkm2_yr => 0.0,
              TROP_Biomass_in_construction_material_being_burnt_GtBiomass_yr => 0.0,
              TROP_Biomass_in_construction_material_left_to_rot_GtBiomass_yr => 0.0,
              TROP_biomass_new_growing_GtBiomass___yr => 0.0,
              TROP_burning_Mkm2_yr => 0.0,
              TROP_Dead_biomass_decomposing_GtBiomass_yr => 0.0,
              TROP_DeadB_and_SOM_tB_per_km2 => 0.0,
              TROP_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr => 0.0,
              TROP_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr => 0.0,
              TROP_deforestation_cutoff => 0.0,
              TROP_deforestation_cutoff_effect => 0.0,
              TROP_deforested_as_pct_of_potial_area => 0.0,
              TROP_deforestion_multiplier_wrt_2000 => 0.0,
              TROP_for_construction_use_GtBiomass_yr => 0.0,
              TROP_historical_deforestation_pct_yr => 0.0,
              TROP_land_taken_out_of_use_GtBiomass => 0.0,
              TROP_land_taken_out_of_use_Mkm2 => 0.0,
              TROP_living_biomass_densitiy_tBiomass_pr_km2 => 0.0,
              TROP_Living_biomass_rotting_GtBiomass_yr => 0.0,
              TROP_NF_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr => 0.0,
              TROP_NF_regrowing_after_being_burnt_Mkm2_yr => 0.0,
              TROP_NF_regrowing_after_harvesting_Mkm2_yr => 0.0,
              TROP_potential_less_actual_living_biomass_GtBiomass => 0.0,
              TROP_potential_living_biomass_GtBiomass => 0.0,
              TROP_regrowing_after_being_clear_cut_Mkm2_yr => 0.0,
              TROP_regrowing_after_being_deforested_Mkm2_yr => 0.0,
              TROP_runoff => 0.0,
              TROP_runoff_time => 0.0,
              TROP_soil_degradation_from_clear_cutting_GtBiomass_yr => 0.0,
              TROP_soil_degradation_from_forest_fires_GtBiomass_yr => 0.0,
              TROP_Sum_outflows_GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass_yr => 0.0,
              TROP_Time_to_decompose_undisturbed_dead_biomass_yr => 0.0,
              TROP_Use_of_NF_biomass_for_energy_GtBiomass_yr => 0.0,
              TROP_with_normal_cover_Mkm2 => 0.0,
              TUNDRA_being_deforested_Mkm2_yr => 0.0,
              TUNDRA_being_harvested_Mkm2_yr => 0.0,
            ]
          end
          push!(startEquationConstructors, generateStartEquations20)
        end
        begin
          function generateStartEquations21()
            [
              TUNDRA_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr => 0.0,
              TUNDRA_Biomass_in_construction_material_being_burnt_GtBiomass_yr => 0.0,
              TUNDRA_Biomass_in_construction_material_left_to_rot_GtBiomass_yr => 0.0,
              TUNDRA_biomass_new_growing_GtBiomass___yr => 0.0,
              TUNDRA_burning_Mkm2_yr => 0.0,
              TUNDRA_Dead_biomass_decomposing_GtBiomass_yr => 0.0,
              TUNDRA_DeadB_and_SOM_tB_per_km2 => 0.0,
              TUNDRA_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr => 0.0,
              TUNDRA_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr => 0.0,
              TUNDRA_for_construction_use_GtBiomass_yr => 0.0,
              TUNDRA_historical_deforestation_pct_yr => 0.0,
              TUNDRA_land_taken_out_of_use_GtBiomass => 0.0,
              TUNDRA_land_taken_out_of_use_Mkm2 => 0.0,
              TUNDRA_living_biomass_densitiy_tBiomass_pr_km2 => 0.0,
              TUNDRA_Living_biomass_rotting_GtBiomass_yr => 0.0,
              TUNDRA_potential_less_actual_living_biomass_GtBiomass => 0.0,
              TUNDRA_potential_living_biomass_GtBiomass => 0.0,
              TUNDRA_regrowing_after_being_burnt_Mkm2_yr => 0.0,
              TUNDRA_regrowing_after_being_deforested_Mkm2_yr => 0.0,
              TUNDRA_regrowing_after_harvesting_Mkm2_yr => 0.0,
              TUNDRA_runoff => 0.0,
              TUNDRA_soil_degradation_from_forest_fires_GtBiomass_yr => 0.0,
              TUNDRA_Sum_outflows_GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass_yr => 0.0,
              TUNDRA_with_normal_cover_Mkm2 => 0.0,
              UNIT_conversion_for_CH4_from_CO2e_to_C => 0.0,
              UNIT_conversion_for_CO2_from_CO2e_to_C => 0.0,
              UNIT_conversion_from_MtCH4_to_GtC => 0.0,
              UNIT_conversion_GtCO2e_to_GtC => 0.0,
              UNIT_conversion_mm_to_m => 0.0,
              UNIT_conversion_W_m2_earth_to_ZJ_yr => 0.0,
              UNIT_converter_GtC_Gm3_to_ymoles_litre => 0.0,
              Upper_to_deep_ocean_temp_diff_in_1850_degC => 0.0,
              Upwelling_from_deep => 0.0,
              Upwelling_to_surface => 0.0,
              Urban_area_fraction => 0.0,
              Urban_Mkm2 => 0.0,
              Urbanzation_Effect_on_biomass_use => 0.0,
              Use_of_GRASS_biomass_for_construction_GtBiomass_yr => 0.0,
              Use_of_GRASS_biomass_for_energy_GtBiomass_yr => 0.0,
              Use_of_GRASS_for_construction_in_2000_GtBiomass => 0.0,
              Use_of_GRASS_for_energy_in_2000_GtBiomass => 0.0,
              Use_of_NF_biomass_for_construction_GtBiomass_yr => 0.0,
              Use_of_NF_biomass_for_energy_GtBiomass_yr => 0.0,
              Use_of_NF_for_construction_in_2000_GtBiomass => 0.0,
              Use_of_NF_for_energy_in_2000_GtBiomass => 0.0,
              Use_of_TROP_biomass_for_construction_GtBiomass_yr => 0.0,
              Use_of_TROP_for_construction_in_2000_GtBiomass => 0.0,
              Use_of_TROP_for_energy_in_2000_GtBiomass => 0.0,
              Use_of_TUNDRA_biomass_for_construction_GtBiomass_yr => 0.0,
              Use_of_TUNDRA_biomass_for_energy_GtBiomass_yr => 0.0,
            ]
          end
          push!(startEquationConstructors, generateStartEquations21)
        end
        begin
          function generateStartEquations22()
            [
              Use_of_TUNDRA_for_construction_in_2000_GtBiomass => 0.0,
              Use_of_TUNDRA_for_energy_in_2000_GtBiomass => 0.0,
              Volcanic_aerosols_emissions => 0.0,
              Volcanic_aerosols_removed_from_stratosphere => 0.0,
              Volume_cold_ocean_0_to_100m => 0.0,
              Volume_cold_ocean_downwelling_100m_to_bottom => 0.0,
              Volume_expansion_from_thermal_expansion_deep_Gm3_km3 => 0.0,
              Volume_expansion_from_thermal_expansion_surface_Gm3_km3 => 0.0,
              Volume_ocean_deep_1km_to_bottom => 0.0,
              Volume_ocean_upwelling_100m_to_1km => 0.0,
              Volume_of_total_ocean_Gm3 => 0.0,
              Volume_warm_ocean_0_to_100m => 0.0,
              Warming_due_to_CH4_blocking_W_m2 => 0.0,
              Warming_due_to_CO2_blocking_W_m2 => 0.0,
              Warming_due_to_othGHG_blocking_W_m2 => 0.0,
              Warming_due_to_water_vapor_blocking_W_m2 => 0.0,
              Years_of_exponential_rise_dless => 0.0,
              Years_of_exponential_rise_yr => 0.0,
              Years_still_needed_to_reach_zero_emission_goal_yr => 0.0,
              yr_on_yr_change_in_C_in_land_use_GtC_yr => 0.0,
              yr_on_yr_change_in_C_in_ocean_GtC_yr => 0.0,
              flow_NF_TUNDRA_Biomass_in_construction_material_left_to_rot_GtBiomass_yr => 0.0,
              flow_NF_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr => 0.0,
              flow_Evaporation_aka_latent_heat_flow => 0.0,
              flow_C_runoff_from_biomass_soil => 0.0,
              flow_Kyoto_Flour_degradation => 0.0,
              flow_N2O_degradation_MtN2O_yr => 0.0,
              flow_LW_TOA_radiation_from_atm_to_space => 0.0,
              flow_TROP_Living_biomass_rotting_GtBiomass_yr => 0.0,
              flow_CO2_flux_TUNDRA_to_atm_Gtc_yr => 0.0,
              flow_Sifting_DESERT_to_GRASS_Mkm2_yr => 0.0,
              flow_Upwelling_from_deep => 0.0,
              flow_TUNDRA_regrowing_after_being_burnt_Mkm2_yr => 0.0,
              flow_TUNDRA_runoff => 0.0,
              flow_Greenland_ice_that_slid_into_the_ocean_melting__pos__or_freezing__neg__km3_yr => 0.0,
              flow_NF_being_harvested_by_clear_cutting_Mkm2_yr => 0.0,
              flow_Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_arctic_ice_ZJ_yr => 0.0,
              flow_TROP_Biomass_in_construction_material_being_burnt_GtBiomass_yr => 0.0,
              flow_Glacial_ice_melting__pos__or_freezing__neg__km3_yr => 0.0,
              flow_NATURE_CCS_Fig3_GtC_yr => 0.0,
              flow_NF_biomass_new_growing_GtBiomass___yr => 0.0,
              flow_LW_clear_sky_emissions_to_surface => 0.0,
              flow_CH4_in_the_atmosphere_converted_to_CO2 => 0.0,
              flow_NF_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr => 0.0,
              flow_TROP_biomass_new_growing_GtBiomass___yr => 0.0,
              flow_GRASS_Living_biomass_rotting_GtBiomass_yr => 0.0,
              flow_TUNDRA_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr => 0.0,
              flow_TUNDRA_biomass_new_growing_GtBiomass___yr => 0.0,
              flow_CO2_flux_from_atm_to_NF_for_new_growth_GtC_yr => 0.0,
              flow_CH4_conversion_to_CO2_and_H2O => 0.0,
            ]
          end
          push!(startEquationConstructors, generateStartEquations22)
        end
        begin
          function generateStartEquations23()
            [
              flow_Flow_of_heat_to_deep_ocean_btw_72_and_08 => 0.0,
              flow_GRASS_for_construction_use_GtBiomass_yr => 0.0,
              flow_Flow_of_cold_surface_water_welling_down_GcubicM_per_yr => 0.0,
              flow_TROP_for_construction_use_GtBiomass_yr => 0.0,
              flow_Annual_glacial_ice_losing__pos__or_gaining__neg__GtIce_yr => 0.0,
              flow_NF_runoff => 0.0,
              flow_NF_soil_degradation_from_forest_fires_GtBiomass_yr => 0.0,
              flow_GRASS_runoff => 0.0,
              flow_Greenland_ice_sliding_into_the_ocean_km3_yr => 0.0,
              flow_TUNDRA_soil_degradation_from_forest_fires_GtBiomass_yr => 0.0,
              flow_Greenland_ice_melting__pos__or_freezing__neg__km3_yr => 0.0,
              flow_SW_surface_absorption => 0.0,
              flow_All_N2O_emissions_MtN2O_yr => 0.0,
              flow_NF_being_harvested_normally_Mkm2_yr => 0.0,
              flow_Kyoto_Flour_emissions => 0.0,
              flow_CO2_flux_from_atm_to_TUNDRA_for_new_growth_GtC_yr => 0.0,
              flow_Shifting_NF_to_TROP_Mkm2_yr => 0.0,
              flow_GRASS_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr => 0.0,
              flow_TUNDRA_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr => 0.0,
              flow_Carbon_flow_from_warm_to_cold_surface_GtC_per_yr => 0.0,
              flow_Shifting_GRASS_to_DESERT_Mkm2_yr => 0.0,
              flow_NF_being_deforested_Mkm2_yr => 0.0,
              flow_Heat_withdrawn_from_atm_by_melting__pos____added__neg__by_freezing_ice_ZJ_yr => 0.0,
              flow_GRASS_biomass_new_growing_GtBiomass___yr => 0.0,
              flow_Man_made_fossil_C_emissions_GtC_yr => 0.0,
              flow_Greenland_ice_melting_as_water_km3_yr => 0.0,
              flow_TROP_runoff => 0.0,
              flow_Antarctic_ice_losing__pos__or_gaining__neg__GtIce_yr => 0.0,
              flow_Carbon_flow_from_cold_surface_downwelling_Gtc_per_yr => 0.0,
              flow_NF_regrowing_after_harvesting_Mkm2_yr => 0.0,
              flow_TROP_Dead_biomass_decomposing_GtBiomass_yr => 0.0,
              flow_TUNDRA_being_deforested_Mkm2_yr => 0.0,
              flow_Shifting_TROP_to_GRASS_Mkm2_yr => 0.0,
              flow_CH4_release___capture_from_permafrost_area_loss___gain_GtC_yr => 0.0,
              flow_Volcanic_aerosols_emissions => 0.0,
              flow_GRASS_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr => 0.0,
              flow_Natural_CH4_emissions => 0.0,
              flow_Flow_of_heat_to_atm_ZJ_yr => 0.0,
              flow_NF_Biomass_in_construction_material_being_burnt_GtBiomass_yr => 0.0,
              flow_Flow_of_heat_to_deep_ocean => 0.0,
              flow_LW_surface_emission => 0.0,
              flow_NF_regrowing_after_being_burnt_Mkm2_yr => 0.0,
              flow_TROP_NF_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr => 0.0,
              flow_TUNDRA_Biomass_in_construction_material_being_burnt_GtBiomass_yr => 0.0,
              flow_Man_made_fossil_C_emissions_for_cumulation_GtC_yr => 0.0,
              flow_C_absorption_by_ocean_from_atm_for_accumulation => 0.0,
              flow_Antarctic_ice_melting__pos__or_freezing__neg__km3_yr => 0.0,
              flow_Annual_flux_of_C_to_biomass_GtC_pr_yr => 0.0,
              flow_GRASS_Biomass_in_construction_material_being_burnt_GtBiomass_yr => 0.0,
              flow_NF_regrowing_after_being_deforested_Mkm2_yr => 0.0,
            ]
          end
          push!(startEquationConstructors, generateStartEquations23)
        end
        begin
          function generateStartEquations24()
            [
              flow_Greenland_ice_losing__pos__or_gaining__neg__GtIce_yr => 0.0,
              flow_CO2_flux_from_atm_to_TROP_for_new_growth_GtC_yr => 0.0,
              flow_NF_soil_degradation_from_clear_cutting_GtBiomass_yr => 0.0,
              flow_Annual_release_of_C_from_permafrost_GtC_y => 0.0,
              flow_Avg_volcanic_activity_GtC_yr => 0.0,
              flow_TUNDRA_regrowing_after_harvesting_Mkm2_yr => 0.0,
              flow_Shifting_ice_on_land_to_tundra_Mkm2_yr => 0.0,
              flow_C_diffusion_into_ocean_from_atm => 0.0,
              flow_Glacial_ice_melting_as_water_km3_yr => 0.0,
              flow_NF_for_construction_use_GtBiomass_yr => 0.0,
              flow_Flow_of_heat_to_surface_ocean => 0.0,
              flow_TROP_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr => 0.0,
              flow_C_released_from_permafrost_as_either_CH4_or_CO2_in_GtC_yr => 0.0,
              flow_TROP_soil_degradation_from_forest_fires_GtBiomass_yr => 0.0,
              flow_TROP_being_harvested_by_clear_cutting_Mkm2_yr => 0.0,
              flow_NF_regrowing_after_being_clear_cut_Mkm2_yr => 0.0,
              flow_GRASS_being_harvested_Mkm2_yr => 0.0,
              flow_Convection_aka_sensible_heat_flow => 0.0,
              flow_TUNDRA_for_construction_use_GtBiomass_yr => 0.0,
              flow_NF_burning_Mkm2_yr => 0.0,
              flow_NF_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr => 0.0,
              flow_TUNDRA_burning_Mkm2_yr => 0.0,
              flow_CO2_flux_TROP_to_atm_GtC_yr => 0.0,
              flow_Shifting_tundra_to_ice_on_land_Mkm2_yr => 0.0,
              flow_Flow_of_water_from_warm_to_cold_surface_Gcubicm_per_yr => 0.0,
              flow_Shifting_Tundra_to_NF_Mkm2_yr => 0.0,
              flow_Flow_of_heat_to_surface_ocean_btw_1972_and_2008 => 0.0,
              flow_TUNDRA_Living_biomass_rotting_GtBiomass_yr => 0.0,
              flow_Methanehydrate_experimental_release_GtC__yr => 0.0,
              flow_GRASS_regrowing_after_being_burnt_Mkm2_yr => 0.0,
              flow_Montreal_gases_degradation => 0.0,
              flow_Carbon_flow_from_cold_to_deep_GtC_per_yr => 0.0,
              flow_GRASS_soil_degradation_from_forest_fires_GtBiomass_yr => 0.0,
              flow_Shifting_TROP_to_NF_Mkm2_yr => 0.0,
              flow_GRASS_being_deforested_Mkm2_yr => 0.0,
              flow_Shifting_GRASS_to_NF_Mkm2_yr => 0.0,
              flow_TROP_being_deforested_Mkm2_yr => 0.0,
              flow_Arctic_ice_melting__pos__or_freezing__neg__km2_yr => 0.0,
              flow_CO2_flux_from_atm_to_GRASS_for_new_growth_GtC_yr => 0.0,
              flow_GRASS_regrowing_after_being_deforested_Mkm2_yr => 0.0,
              flow_Net_C_to_atm_rate => 0.0,
              flow_Methane_hydrates_released_and_converted_to_CO2_by_bacteria_GtC => 0.0,
              flow_LW_surface_emissions_NOT_escaping_through_atm_window => 0.0,
              flow_Antarctic_ice_melting_as_water_km3_yr => 0.0,
              flow_TUNDRA_Biomass_in_construction_material_left_to_rot_GtBiomass_yr => 0.0,
              flow_TROP_NF_regrowing_after_harvesting_Mkm2_yr => 0.0,
              flow_TUNDRA_being_harvested_Mkm2_yr => 0.0,
              flow_Net_heat_flow_ocean_from_surface_to_deep__ZJ_yr_ => 0.0,
              flow_TROP_regrowing_after_being_clear_cut_Mkm2_yr => 0.0,
              flow_GRASS_Biomass_in_construction_material_left_to_rot_GtBiomass_yr => 0.0,
            ]
          end
          push!(startEquationConstructors, generateStartEquations24)
        end
        begin
          function generateStartEquations25()
            [
              flow_TROP_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr => 0.0,
              flow_Carbon_flow_from_deep => 0.0,
              flow_Rate_of_destruction_of_wetlands => 0.0,
              flow_Montreal_gases_emissions => 0.0,
              flow_LW_re_radiated_by_clouds => 0.0,
              flow_Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_antarctic_ice_ZJ_yr => 0.0,
              flow_Depositing_of_C_to_sediment => 0.0,
              flow_TUNDRA_Dead_biomass_decomposing_GtBiomass_yr => 0.0,
              flow_TUNDRA_regrowing_after_being_deforested_Mkm2_yr => 0.0,
              flow_TROP_burning_Mkm2_yr => 0.0,
              flow_TROP_NF_regrowing_after_being_burnt_Mkm2_yr => 0.0,
              flow_SW_Atmospheric_absorption => 0.0,
              flow_GRASS_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr => 0.0,
              flow_GRASS_regrowing_after_harvesting_Mkm2_yr => 0.0,
              flow_TROP_being_harvested_normally_Mkm2_yr => 0.0,
              flow_C_release_from_permafrost_melting_as_CO2_GtC_yr => 0.0,
              flow_Human_activity_CH4_emissions => 0.0,
              flow_GRASS_Dead_biomass_decomposing_GtBiomass_yr => 0.0,
              flow_TROP_Biomass_in_construction_material_left_to_rot_GtBiomass_yr => 0.0,
              flow_TROP_soil_degradation_from_clear_cutting_GtBiomass_yr => 0.0,
              flow_TUNDRA_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr => 0.0,
              flow_Shifting_NF_to_GRASS_Mkm2_yr => 0.0,
              flow_Heat_flow_from_the_earths_core => 0.0,
              flow_Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_Greenland_ice_that_slid_into_the_ocean_ZJ_yr => 0.0,
              flow_TROP_regrowing_after_being_deforested_Mkm2_yr => 0.0,
              flow_C_removal_rate_from_atm_for_nature_May_2020_GtC_y => 0.0,
              flow_GRASS_burning_Mkm2_yr => 0.0,
              flow_CO2_flux_GRASS_to_atm_Gtc_yr => 0.0,
              flow_Upwelling_to_surface => 0.0,
              flow_NF_Dead_biomass_decomposing_GtBiomass_yr => 0.0,
              flow_Carbon_captured_and_stored_GtC___yr => 0.0,
              flow_Volcanic_aerosols_removed_from_stratosphere => 0.0,
              flow_Carbon_flow_from_intermediate_to_surface_box_GtC_per_yr => 0.0,
              flow_Greenland_ice_melting_that_slid_into_the_ocean_km3_yr => 0.0,
              flow_Shifting_NF_to_Tundra_Mkm2_yr => 0.0,
              flow_Shifting_GRASS_to_TROP_Mkm2_yr => 0.0,
              flow_NF_Living_biomass_rotting_GtBiomass_yr => 0.0,
              flow_CO2_flux_NF_to_atm_Gtc_yr => 0.0,
              flow_Flow_of_cold_water_sinking_to_very_bottom_GcubicM_per_yr => 0.0,
              flow_Biological_removal_of_C_from_WSW_GtC_per_yr => 0.0,
              C_in_ocean_1_yr_ago_GtC_DL => 0.0,
              Atmos_heat_used_for_melting_last_year_1_yr => 0.0,
              Ocean_heat_used_for_melting_last_year_ZJ_yr => 0.0,
              C_in_atm_1_yr_ago_GtC => 0.0,
              C_in_atm_1_yr_ago_GtC_RT1 => 0.0,
              C_in_atm_1_yr_ago_GtC_RT2 => 0.0,
              C_in_atm_1_yr_ago_GtC_DL => 0.0,
              All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC => 0.0,
              All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_RT1 => 0.0,
              All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_RT2 => 0.0,
            ]
          end
          push!(startEquationConstructors, generateStartEquations25)
        end
        begin
          function generateStartEquations26()
            [
              All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_DL => 0.0,
              ifEq_tmp304 => 0.0,
              ifEq_tmp305 => 0.0,
              Arctic_land_surface_temp_anomaly_compared_to_1850 => 0.0,
              Biological_removal_of_C_from_WSW_GtC_per_yr => 0.0,
              Effect_of_temp_on_permafrost_melting_dmnl => 0.0,
              Temp_diff_relevant_for_melting_or_freezing_arctic_ice_from_1850 => 0.0,
              Temp_diff_relevant_for_melting_or_freezing_from_1850 => 0.0,
              yr_on_yr_change_in_C_in_atm_GtC_yr => 0.0,
              C_in_ocean_1_yr_ago_GtC => 0.0,
              C_in_ocean_1_yr_ago_GtC_LV1 => 0.0,
              C_in_ocean_1_yr_ago_GtC_LV2 => 0.0,
              Atmos_heat_used_for_melting_last_year_1_yr_LV => 0.0,
              Ocean_heat_used_for_melting_last_year_ZJ_yr_LV => 0.0,
              C_in_atm_1_yr_ago_GtC_LV1 => 0.0,
              C_in_atm_1_yr_ago_GtC_LV2 => 0.0,
              C_in_atm_1_yr_ago_GtC_LV3 => 0.0,
              All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_LV1 => 0.0,
              All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_LV2 => 0.0,
              All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_LV3 => 0.0,
              Cumulative_carbon_removed_from_atm_for_nature_May_2020 => 0.0,
              Antarctic_ice_volume_km3 => 3.0e7,
              Arctic_ice__on_sea__area_km2 => 1.34e7,
              C_in_atmosphere_GtC => 600.0,
              C_in_atmosphere_in_form_of_CH4 => 1.69,
              C_in_cold_surface_water_GtC => Carbon_in_cold_ocean_0_to_100m_1850_GtC,
              C_in_cold_water_trunk_downwelling_GtC => Carbon_in_cold_ocean_trunk_100m_to_bottom_1850_GtC,
              C_in_deep_water_volume_1km_to_bottom_GtC => Carbon_in_ocean_deep_1k_to_bottom_ocean_1850_GtC,
              C_in_intermediate_upwelling_water_100m_to_1km_GtC => Carbon_in_ocean_upwelling_100m_to_1km_1850_GtC,
              C_in_permafrost_in_form_of_CH4 => 1200.0,
              C_in_sediment => 3.0e9,
              C_in_warm_surface_water_GtC => Carbon_in_warm_ocean_0_to_100m_1850_GtC,
              Cold_surface_water_volume_Gm3 => Volume_cold_ocean_0_to_100m,
              Cold_water_volume_downwelling_Gm3 => Volume_cold_ocean_downwelling_100m_to_bottom,
              Cumulative_antarctic_ice_volume_loss_GtIce => 0.0,
              Cumulative_C_released_from_permafrost_as_either_CH4_or_CO2_in_GtC_yr => 0.0,
              Cumulative_carbon_captured_and_stored_GtC => 0.0,
              Cumulative_carbon_removed_from_atm_for_nature_May_2020 => 0.0,
              Cumulative_flow_of_C_to_biomass_since_1850_GtC => 0.0,
              Cumulative_glacial_ice_volume_loss_GtIce => 0.0,
              Cumulative_Greenland_ice_volume_loss_GtIce => 0.0,
              Cumulative_heat_to_atm_ZJ => 0.0,
              Cumulative_ocean_volume_increase_due_to_ice_melting_km3 => 0.0,
              Cumulative_release_of_C_from_permafrost_GtC => 0.0,
              Deep_water_volume_1km_to_4km_Gm3 => Volume_ocean_deep_1km_to_bottom,
              DESERT_Mkm2 => 25.4,
              Fossil_fuel_reserves_in_ground_GtC => 6000.0,
              Glacial_ice_volume_km3 => 167000.0,
              GRASS_area_burnt_Mkm2 => 1.0,
              GRASS_area_harvested_Mkm2 => 2.5,
            ]
          end
          push!(startEquationConstructors, generateStartEquations26)
        end
        begin
          function generateStartEquations27()
            [
              GRASS_Biomass_locked_in_construction_material_GtBiomass => 1.5,
              GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass => 1200.0,
              GRASS_deforested_Mkm2 => 0.5,
              GRASS_Living_biomass_GtBiomass => 310.0,
              GRASS_potential_area_Mkm2 => 22.5,
              Greenland_ice_volume_on_Greenland_km3 => 2.93e6,
              Greenland_ice_volume_that_slid_into_the_ocean_km3 => 0.0,
              Heat_in_atmosphere_ZJ => 1025.67,
              Heat_in_deep_ZJ => 1.9532e6,
              Heat_in_surface => 25000.0,
              Intermediate_upwelling_water_volume_100m_to_1km_Gm3 => Volume_ocean_upwelling_100m_to_1km,
              Kyoto_Flour_gases_in_atm => 0.0,
              Montreal_gases_in_atm => 0.0,
              N2O_in_atmosphere_MtN2O => 900.0,
              NATURE_Cumulative_CCS_GtC => 0.0,
              NF_area_burnt_Mkm2 => 2.5,
              NF_area_clear_cut_Mkm2 => 1.0,
              NF_area_deforested_Mkm2 => 0.0,
              NF_area_harvested_Mkm2 => 1.0,
              NF_Biomass_locked_in_construction_material_GtBiomass => 3.0,
              NF_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass => 330.0,
              NF_Living_biomass_GtBiomass => 115.0,
              NF_potential_area_Mkm2 => 17.0,
              Sum_C_absorbed_by_ocean_GtC => 0.0,
              Sum_heat_to_deep_ocean => 0.0,
              Sum_heat_to_deep_ocean_btw_72_and_08 => 0.0,
              Sum_heat_to_surface_ocean_btw_72_and_08 => 0.0,
              Sum_heat_to_surface_ocean_ZJ => 0.0,
              Sum_man_made_CO2_emissions_GtC => 0.0,
              Sum_net_C_to_atm => 0.0,
              TROP_area_burnt_Mkm2 => 1.7,
              TROP_area_clear_cut_Mkm2 => 0.3,
              TROP_area_deforested_Mkm2 => 1.0,
              TROP_area_harvested_Mkm2 => 0.3,
              TROP_Biomass_locked_in_construction_material_GtBiomass => 30.0,
              TROP_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass => 160.0,
              TROP_Living_biomass_GtBiomass => 370.0,
              TROP_potential_area_Mkm2 => 25.0,
              TUNDRA_area_burnt_Mkm2 => 2.0,
              TUNDRA_area_harvested_Mkm2 => 2.5,
              TUNDRA_Biomass_locked_in_construction_material_GtBiomass => 1.5,
              TUNDRA_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass => 1200.0,
              TUNDRA_deforested_Mkm2 => 0.0,
              TUNDRA_Living_biomass_GtBiomass => 300.0,
              Tundra_potential_area_Mkm2 => 22.5,
              Volcanic_aerosols_in_stratosphere => 0.0,
              Warm_surface_water_volume_Gm3 => Volume_warm_ocean_0_to_100m,
              Wetlands_area => 1.0e7,
              Aerosol_anthropogenic_emissions_in_2010 => 0.0,
              CO2_emissions_in_2010 => 0.0,
            ]
          end
          push!(startEquationConstructors, generateStartEquations27)
        end
        begin
          function generateStartEquations28()
            [
              CO2_ppm_value_at_When_to_sample => MODEL_CO2_concentration_in_atmosphere2_ppm,
              CO4_emissions_in_2010 => 0.0,
              Greenland_slide_experiment_end_condition => 0.0,
              Kyoto_Flour_concentration_in_1970_ppt => 0.0,
              Kyoto_Flour_emissions_RCPs_JR_in_2010 => 0.0,
              Montreal_gases_concentration_in_1970_ppt => 0.0,
              Montreal_gases_emissions_RCPs_JR_in_2010 => 0.0,
              N20_emissions_RCPs_JR_in_2010 => 0.0,
              Tipping_point_search_amount_at_start => 12.0,
              Arctic_land_surface_temp_anomaly_compared_to_1850 => Temp_surface_anomaly_compared_to_1850_degC,
              Biological_removal_of_C_from_WSW_GtC_per_yr => Net_marine_primary_production_NMPP_GtC_pr_yr,
              Effect_of_temp_on_permafrost_melting_dmnl => 1.0 + Slope_btw_temp_and_permafrost_melting___freezing * (Temp_diff_relevant_for_melting_or_freezing_from_1850 / 4.0 - 1.0),
              Temp_diff_relevant_for_melting_or_freezing_arctic_ice_from_1850 => Temp_surface_anomaly_compared_to_1850_degC,
              Temp_diff_relevant_for_melting_or_freezing_from_1850 => Temp_surface_C - 13.66500000000002,
              yr_on_yr_change_in_C_in_atm_GtC_yr => C_in_atmosphere_GtC - C_in_atm_1_yr_ago_GtC,
              C_in_ocean_1_yr_ago_GtC => Total_carbon_in_ocean_GtC,
              C_in_ocean_1_yr_ago_GtC_LV1 => Total_carbon_in_ocean_GtC,
              C_in_ocean_1_yr_ago_GtC_LV2 => Total_carbon_in_ocean_GtC,
              Atmos_heat_used_for_melting_last_year_1_yr_LV => 0.0,
              Ocean_heat_used_for_melting_last_year_ZJ_yr_LV => 0.0,
              C_in_atm_1_yr_ago_GtC_LV3 => C_in_atm_1_yr_ago_GtC_DL * C_in_atmosphere_GtC,
              C_in_atm_1_yr_ago_GtC_LV2 => C_in_atm_1_yr_ago_GtC_LV3,
              C_in_atm_1_yr_ago_GtC_LV1 => C_in_atm_1_yr_ago_GtC_LV3,
              All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_LV3 => All_C_taken_out_due_to_change_in_land_use_GtC * All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_DL,
              All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_LV2 => All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_LV3,
              All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_LV1 => All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_LV3,
            ]
          end
          push!(startEquationConstructors, generateStartEquations28)
        end
      end
      for constructor in startEquationConstructors
        push!(startEquationComponents, constructor())
      end
      initialValues = collect(Iterators.flatten(startEquationComponents))
      equationComponents = []
      begin
        begin
          Future_volcanic_emissions = 0.0
          Albedo_Antarctic_hist = 0.7
          Albedo_Antarctic_sens = 0.7
          Albedo_BARREN_normal = 0.17
          Albedo_BARREN_white = 0.7
          Albedo_DESERT_normal = 0.24
          Albedo_glacier_hist = 0.4
          Albedo_glacier_sens = 0.4
          Albedo_GRASS_burnt = 0.08
          Albedo_GRASS_deforested = 0.3
          Albedo_GRASS_normal_cover = 0.16
          Albedo_Greenland = 0.7
          Albedo_NF_burnt = 0.13
          Albedo_NF_deforested = 0.18
          Albedo_NF_normal_cover = 0.08
          Albedo_TROP_burnt = 0.1
          Albedo_TROP_deforested = 0.168
          Albedo_TROP_normal_cover = 0.14
          Albedo_TUNDRA_burnt = 0.23
          Albedo_TUNDRA_deforested = 0.23
          Albedo_TUNDRA_normal_cover = 0.23
          Albedo_URBAN_normal = 0.15
          Amount_methane_hydrates__clathrates__experimentally_released_GtC = 0.0
          Amt_of_constant_emissions_GtC_yr = 4.0
          Annual_pct_increase_CH4_emissions_from_2015_pct_yr = 0.0
          Annual_pct_increase_CO2_emissions_from_2015_pct_yr = 0.0
          Antarctic_ice_volume_in_1850_km3 = 3.0e7
          Arctic_ice_albedo_1850 = 0.7
          Arctic_ice_area_in_1850_km2 = 1.34e7
          Arctic_surface_temp_delay_yr = 15.0
          Area_covered_by_high_clouds_in_1850 = 0.2
          Area_covered_by_low_clouds_in_1850 = 0.4
          Area_equivalent_of_1km_linear_retreat_km2 = 17500.0
          Area_of_earth_m2 = 5.1e14
          Area_of_ocean_at_surface_361900_Gm2 = 361900.0
          Atmos_heat_used_for_melting_Initially_1_yr = 0.0
          Average_thickness_arctic_ice_km = 0.0025
          Avg_amount_of_C_in_the_form_of_CH4_per_km2 = 4.8e-5
          Avg_depth_of_permafrost_km = 0.1
          Avg_flatness_of_worlds_coastline = 1.0
          Avg_thickness_Antarctic_hist_km = 2.14
          Avg_thickness_Antarctic_sens_km = 2.14
          Avg_thickness_Greenland_km = 1.35
          C_in_atmosphere_in_1850_GtC = 600.0
          C_in_the_form_of_CH4_in_atm_1850 = 1.69
          Carbon_per_biomass_tC_per_tBiomass = 0.5
          CC_in_cold_ocean_0_to_100m_1850_ymoles_per_litre = 2240.0
          CC_in_cold_ocean_downwelling_100m_bottom_1850_ymoles_per_litre = 2240.0
          CC_in_ocean_upwelling_100m_to_1km_1850_ymoles_per_litre = 2240.0
          CC_in_warm_ocean_0_to_100m_1850_ymoles_per_litre = 2240.0
          CC_ocean_deep_1km_to_bottom_1850_ymoles_per_litre = 2240.0
          CH4_concentration_in_2010_ppb = 1720.81
          CH4_halflife_in_atmosphere = 7.3
          Cold_dense_water_sinking_in_Sverdrup_in_1850 = 35.0
          Constant_anthropogenic_CH4_emissions = 0.2
          Convection_as_f_of_incoming_solar_in_1850 = 0.071
          conversion_factor_CH4_Gt_to_ppb = 468.0
          Conversion_from_Kyoto_Flour_amount_to_concentration_ppt_kt = 0.04
          Conversion_from_Montreal_gases_amount_to_concentration_ppt_kt = 0.04
          Conversion_Millionkm2_to_km2_Mkm2_km2 = 1.0e-6
          Conversion_of_anthro_aerosol_emissions_to_forcing = -1.325
          Conversion_of_volcanic_aerosol_emissions_to_CO2_emissions_GtC_pr_VAE = 2.8
          Conversion_of_volcanic_aerosol_forcing_to_volcanic_aerosol_emissions = -1.0
          Conversion_ymoles_per_kg_to_pCO2_yatm = 0.127044
          Densitiy_of_water_relative_to_ice = 0.916
          Duration_of_destruction_yr = 5.0
          Emissions_of_natural_CH4_GtC_yr = 0.19
          Emissivity_atm = 1.0
          Emissivity_surface = 1.0
          Evaporation_as_fraction_of_incoming_solar_in_1850 = 0.289
          EXP_12f_Stratospheric_scattering_experiment_0_off_1_on = float(0)
          Experimental_doubling_of_constant_C_emissions_how_long_yr = 5.0
          Experimental_doubling_of_constant_C_emissions_how_much_1_100pct = 0.0
          Experimental_doubling_of_constant_C_emissions_when_yr = 30000.0
          Frac_of_surface_emission_through_atm_window = 0.051
          Frac_SW_clear_sky_reflection_aka_scattering = 0.0837
          Frac_SW_HI_cloud_efffect_aka_cloud_albedo = 0.006
          Frac_SW_LO_cloud_efffect_aka_cloud_albedo = 0.158
          Fraction_of_C_released_from_permafrost_released_as_CH4_hist_dmnl = 1.0
          Fraction_of_C_released_from_permafrost_released_as_CH4_sensitivity_dmnl = 1.0
          Fraction_of_earth_surface_as_ocean = 0.7
          Fraction_of_heat_needed_to_melt_antarctic_ice_coming_from_air = 0.6
          Fraction_of_heat_needed_to_melt_arctic_ice_coming_from_air = 0.5
          Fraction_of_heat_needed_to_melt_Greenland_ice_that_slid_into_the_ocean_coming_from_air = 0.1
          Fraction_of_methane_hydrates_released_from_the_ocean_converted_to_CO2_before_it_is_relased_to_the_atmosphere = 0.9
          Fraction_of_ocean_classified_warm_surface = 0.8
          Glacial_ice_volume_in_1850_km3 = 167000.0
          Global_Warming_Potential_CH4 = 25.0
          Global_Warming_Potential_N20 = 298.0
          GRASS_area_burned_in_1850_Mkm2 = 1.0
          GRASS_area_deforested_in_1850_Mkm2 = 0.5
          GRASS_area_harvested_in_1850_Mkm2 = 2.5
          GRASS_Avg_life_biomass_yr = 100.0
          GRASS_Avg_life_of_building_yr = 10.0
          GRASS_Biomass_locked_in_construction_material_in_1850_GtBiomass = 1.5
          GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_in_1850_GtBiomass = 1200.0
          GRASS_Fraction_of_construction_waste_burned_0_1 = 0.5
          GRASS_fraction_of_DeadB_and_SOM_being_destroyed_by_deforesting = 1.0
          GRASS_fraction_of_DeadB_and_SOM_being_destroyed_by_energy_harvesting = 0.1
          GRASS_fraction_of_DeadB_and_SOM_destroyed_by_natural_fires = 0.0
          GRASS_living_biomass_densitiy_in_1850_tBiomass_pr_km2 = 14500.0
          GRASS_Living_biomass_in_1850_GtBiomass = 310.0
          GRASS_Normal_fire_incidence_fraction_yr = 1.0
          GRASS_Ref_historical_deforestation_pct_yr = 0.1
          GRASS_runoff_time = 2000.0
          GRASS_Speed_of_regrowth_yr = 2.0
          GRASS_Time_to_decompose_undisturbed_dead_biomass_yr = 1000.0
          Greenland_ice_slide_circulation_slowdown_effect = 0.33
          Greenland_ice_volume_in_1850_km3 = 2.93e6
          Greenland_slide_experiment_how_much_sildes_in_the_ocean_fraction = 0.25
          Greenland_slide_experiment_over_how_many_years_yr = 70.0
          GtIce_vs_km3 = 0.9167
          Heat_gained___needed_to_freeze___unfreeze_1_km3_permafrost_ZJ_km3 = 0.0001717
          Heat_in__ocean__deep_in_1850_ZJ = 1.9532e6
          Heat_in_atmosphere_in_1850_ZJ = 1025.67
          Heat_in_surface_in_1850_ZJ = 25000.0
          Heat_needed_to_melt_1_km3_of_ice_ZJ = 0.0003327
          Hist_Avg_thickness_glacier_km = 0.23
          Hist_Net_heat_flow_ocean_between_surface_and_deep_per_K_of_difference_ZJ_yr_K = 10.0
          Hist_NF_Avg_life_biomass_yr = 60.0
          Hist_NF_Speed_of_regrowth_yr = 3.0
          Hist_Slope_of_effect_of_temp_shifting_DESERT_to_GRASS = 0.4
          Hist_Slope_temp_vs_glacial_ice_melting = 1.0
          Hist_Time_in_trunk = 234.638
          Hist_Time_to_degrade_Kyoto_Flour_yr = 50.0
          Hist_Time_to_regrow_NF_after_buning_yr = 30.0
          Hist_TROP_runoff_time = 2000.0
          Hist_TROP_Time_to_decompose_undisturbed_dead_biomass_yr = 24.0
          K_to_C_conversion_C_K = 273.15
          Kyoto_Flour_Global_Warming_Potential = 7000.0
          Land_surface_temp_adjustment_time_yr = 25.0
          LW_ALL_cloud_radiation_reference_in_1850_W_m2 = 27.9
          LW_LO_cloud_radiation_reference_in_1850_W_m2 = 20.0
          LW_radiation_fraction_blocked_by_other_GHG_in_1850 = 0.0398
          Man_made_CH4_emissions_in_2015_GtC = 0.303
          Man_made_CO2_emissions_in_2015_GtC = 10.0
          MAX_NATURE_CCS_removal_in_2050_GtCO2e_yr = 35.0
          Melting_of_permafrost_at_all_depths_at_4_deg_C_temp_diff_km_yr = 0.71
          Montreal_Global_Warming_Potential = 10000.0
          Myhre_constant_for_CH4 = 0.0594
          Myhre_constant_for_CO2 = 5.35
          Myhre_constant_for_N20 = 0.12
          N2O_concentration_in_2010_ppb = 363.504
          N2O_in_atmosphere_MtN2O_in_1850 = 900.0
          N2O_natural_emissions = 9.0
          Net_marine_primary_production_in_1850 = 0.4
          NEvt_13a_double_rate_of_melting_ice_and_permafrost = float(1)
          NEvt_13b2_Double_incidence_of_biomass_fires = float(1)
          NEvt_13b3_double_sunspot_amplitude_from_2015_onwards_1_normal_2_double = float(1)
          NEvt_13c1_increase_in_area_covered_by_low_clouds = float(1)
          NEvt_13d_Greenland_slide_experiment_start_yr = float(3000000)
          NEvt_2a_Volcanic_eruptions_in_the_future_VAEs_first_future_pulse = float(21000000)
          NEvt_3b_increase_in_area_covered_by_high_clouds = float(1)
          NF_area_burned_in_1850_Mkm2 = 2.5
          NF_area_deforested_in_1850_Mkm2 = 0.0
          NF_area_harvested_in_1850_Mkm2 = 1.0
          NF_Avg_life_of_building_yr = 20.0
          NF_Biomass_locked_in_construction_material_in_1850_GtBiomass = 3.0
          NF_Dead_biomass__litter_and_soil_organic_matter_SOM_in_1850_GtBiomass = 330.0
          NF_DeadB_and_SOM_densitiy_in_1850_tBiomass_pr_km2 = 27500.0
          NF_Fraction_of_construction_waste_burned_0_1 = 0.5
          NF_fraction_of_DeadB_and_SOM_being_destroyed_by_clear_cutting = 0.5
          NF_fraction_of_DeadB_and_SOM_being_destroyed_by_deforesting = 1.0
          NF_fraction_of_DeadB_and_SOM_being_destroyed_by_energy_harvesting = 0.1
          NF_fraction_of_DeadB_and_SOM_destroyed_by_natural_fires = 0.0
          NF_living_biomass_densitiy_in_1850_tBiomass_pr_km2 = 7500.0
          NF_Living_biomass_in_1850_GtBiomass = 115.0
          NF_Normal_fire_incidence_fraction_yr = 0.7
          NF_Ref_historical_deforestation_pct_yr = 0.02
          NF_runoff_time = 2000.0
          NF_Time_to_decompose_undisturbed_dead_biomass_yr = 250.0
          Ocean_heat_used_for_melting_Initially_1_yr = 0.0
          Ocean_slowdown_experimental_factor = 1.0
          Open_ocean_albedo = 0.065
          Over_how_many_yrs_methane_hydrate_release_yr = 5.0
          per_annum_yr = 1.0
          Policy_1_Reducing_GHG_emissions_by_one_third_by_2035 = float(0)
          Policy_2_Large_scale_implementation_of_carbon_capture_and_geological_storage__CCS_ = float(0)
          Population_2000_bn = 6.1
          Pressure_adjustment_deep_pct = 1.0
          Pressure_adjustment_surface_pct = 0.2
          Rate_of_wetland_destruction_pct_of_existing_wetlands_yr = 0.0
          Ratio_of_methane_in_tundra_to_wetland = 4.0
          Ref_shifting_biome_yr = 50.0
          Ref_temp_difference__4_degC_ = 4.0
          Ref_temp_difference_for_antarctic_ice_melting__3_degC_ = 3.0
          Ref_temp_difference_for_Arctic_ice_melting = 0.4
          Ref_temp_difference_for_glacial_ice_melting__1_degC_ = 3.0
          Ref_temp_difference_for_greenland_ice_melting_C = 1.0
          Ref_temp_difference_for_greenland_ice_that_slid_into_the_ocean_melting_degC = 1.0
          Reference_temp_C = 10.0
          Reference_Time_to_regrow_TROP_after_deforesting_yr = 10000.0
          SCALE_and_UNIT_converter_zero_C_to_K = 273.15
          Sens_Avg_thickness_glacier_km = 0.23
          Sens_Frac_atm_absorption = 0.220588
          Sens_Net_heat_flow_ocean_between_surface_and_deep_per_K_of_difference_ZJ_yr_K = 10.0
          Sens_NF_Avg_life_biomass_yr = 60.0
          Sens_NF_Speed_of_regrowth_yr = 3.0
          Sens_Slope_of_effect_of_temp_shifting_DESERT_to_GRASS = 0.4
          Sens_Slope_temp_vs_glacial_ice_melting = 1.0
          Sens_Time_in_trunk = 234.638
          Sens_Time_to_degrade_Kyoto_Flour_yr = 50.0
          Sens_Time_to_regrow_NF_after_buning_yr = 30.0
          Sens_TROP_runoff_time = 2000.0
          Sens_TROP_Time_to_decompose_undisturbed_dead_biomass_yr = 24.0
          Sensitivity_of_biomass_new_growth_to_CO2_concentration = 1.0
          Sensitivity_of_convection_to_temp = 2.5
          Sensitivity_of_evaporation_to_temp = 0.58
          Sensitivity_of_high_cloud_coverage_to_temp_base = 50.0
          Sensitivity_of_high_cloud_coverage_to_temp_sens = 50.0
          Sensitivity_of_low_cloud_coverage_to_temp = 58.0
          Sensitivity_of_trop_to_humidity = 5.0
          Slider_for_annual_removal_of_C_from_atm_after_2020_GtC_y = 0.0
          Slider_for_H2O_slope_hist = 0.0
          Slider_for_slope_fut = 0.0
          Slope_btw_Kyoto_Flour_ppt_and_blocking_multiplier = 0.3
          Slope_btw_Montreal_gases_ppt_and_blocking_multiplier = 0.3
          Slope_btw_N2O_ppb_and_blocking_multiplier = 0.1
          Slope_btw_temp_and_permafrost_melting___freezing_base = 1.0
          Slope_btw_temp_and_permafrost_melting___freezing_sensitivity = 1.0
          Slope_Effect_Temp_on_NMPP = 2.0
          Slope_of_effect_of_temp_on_shifting_NF_to_Tundra = 0.1
          Slope_of_effect_of_temp_on_shifting_TROP_to_NF = 1.0
          Slope_of_effect_of_temp_shifting_GRASS_to_DESERT = 5.0
          Slope_of_effect_of_temp_shifting_GRASS_to_NF = 0.1
          Slope_of_effect_of_temp_shifting_GRASS_to_TROP = 0.2
          Slope_of_effect_of_temp_shifting_NF_to_GRASS = 0.01
          Slope_of_effect_of_temp_shifting_NF_to_TROP = 0.2
          Slope_of_effect_of_temp_shifting_TROP_to_GRASS = 0.05
          Slope_of_effect_of_temp_shifting_tundra_to_NF = 0.2
          Slope_of_efffect_of_acidification_on_NMPP = 5.0
          Slope_temp_eff_on_fire_incidence = 0.1
          Slope_temp_vs_antarctic_ice_melting = 1.2
          Slope_temp_vs_Arctic_ice_melting = 0.65
          Slope_temp_vs_greenland_ice_melting = 0.1
          Slope_temp_vs_greenland_ice_that_slid_into_the_ocean_melting = 0.71
          Solar_sine_forcing_amplitude = 0.1
          Solar_sine_forcing_lift = 0.05
          Solar_sine_forcing_offset_yr = -3.5
          Solar_sine_forcing_period_yr = 11.0
          Stephan_Boltzmann_constant = 5.67037e-8
          Stratospheric_scattering_experiment_end_year = 3.0e7
          Stratospheric_scattering_experiment_reduction_from_2015_in_W_m2 = 3.0
          Switch_0_normal_model_1_dbl_CO2_2_1pct_incr = float(0)
          Switch_btw_historical_CO2_CH4_emissions_or_constant_1history_0constant = float(1)
          SWITCH_for_NATURE_comm_200115_base_1_cut_all_mm_emi_in_2020_2 = float(1)
          SWITCH_future_slope_base_0_plus_5_1_minus_5_2 = float(0)
          SWITCH_h2o_blocked_table_0_linear_1_poly_2 = float(2)
          SWITCH_h2o_poly_dyn_0_equ_1 = float(1)
          SWITCH_nature_rev_0_base_1_steeper_2_less_steep = float(0)
          Switch_to_choose_CO2_CH4_aerosol_other_GHG_0normal_1zero_from_2010_2constant_from_2010 = float(0)
          Switch_to_choose_input_emission_scenario_for_CO2_CH4_and_oth_GHG = float(1)
          Switch_to_drive_model_with_normal_ESCIMO_data__0__CO2e_from_C_Roads__1__or_CO2e_from_CAT_2__or_user_determined_CO2_max_to_find_temp_tipping_point__3_ = 0.0
          Switch_to_run_experiment_12a_reduction_in_emissions_0_off_1_on = float(0)
          Switch_to_run_experiment_12b_CCS_0_off_1_on = float(0)
          Switch_to_run_experiment_12c_stopping_TROP_deforestation_0_off_1_on = float(0)
          Switch_to_run_experiment_12e_white_surfaces_0_off_1_on = float(0)
          Switch_to_run_NATURE_experiment_CCS_0_off_1_on_0 = float(0)
          Switch_to_run_POLICY_4_Stopping_logging_in_Northern_forests_0_off_1_on = float(0)
          Temp__ocean__deep_in_1850_C = 4.0
          Temp_atm_1850 = 274.31
          Temp_gradient_in_surface_degK = 9.7
          Temp_surface_1850_K = 286.815
          TEST_Year_in_which_zero_emissions_are_to_be_reached_yr_Remember_to_set_switch_to_9Linear = 2050.0
          Thickness_of_deep_water_box_1km_to_bottom = 2800.0
          Thickness_of_intermediate_water_box_800m = 800.0
          Thickness_of_surface_water_box_100m = 100.0
          Time_at_which_human_deforestation_is_stopped = 3000.0
          Time_for_volcanic_aerosols_to_remain_in_the_stratosphere = 1.0
          Time_in_cold = 6.51772
          Time_in_deep = 739.89
          Time_in_intermediate_yr = 211.397
          Time_in_warm = 26.227
          Time_to_degrade_Montreal_gases_yr = 30.0
          Time_to_degrade_N2O_in_atmopshere_yr = 95.0
          Time_to_deposit_C_in_sediment = 20000.0
          Time_to_let_shells_form_and_sink_to_sediment_yr = 25.0
          Time_to_melt_Arctic_ice_at_the_reference_delta_temp = 500.0
          Time_to_melt_greenland_ice_at_the_reference_delta_temp = 4000.0
          Time_to_melt_greenland_ice_that_slid_into_the_ocean_at_the_reference_delta_temp = 20.0
          Time_to_melt_or_freeze_antarctic_ice_at_the_reference_delta_temp = 18000.0
          Time_to_melt_or_freeze_glacial_ice_at_the_reference_delta_temp = 500.0
          Time_to_propagate_temperature_change_through_the_volume_of_permafrost_yr = 5.0
          Time_to_reach_C_equilibrium_between_atmosphere_and_ocean = 18.0
          Time_to_regrow_GRASS_after_buning_yr = 10.0
          Time_to_regrow_GRASS_after_deforesting_yr = 80.0
          Time_to_regrow_NF_after_deforesting_yr = 80.0
          Time_to_regrow_TROP_after_buning_yr = 30.0
          Time_to_regrow_TUNDRA_after_buning_yr = 10.0
          Time_to_regrow_TUNDRA_after_deforesting_yr = 80.0
          Time_to_smooth_out_temperature_diff_relevant_for_melting_or_freezing_from_1850_yr = 3.0
          Tipping_point_search_amount_at_peak = 0.0
          Tipping_point_year_of_end = 210000.0
          Tipping_point_year_of_start = 500000.0
          TROP_area_burned_in_1850_Mkm2 = 1.7
          TROP_area_deforested_in_1850_Mkm2 = 1.0
          TROP_area_harvested_in_1850_Mkm2 = 0.3
          TROP_Avg_life_biomass_yr = 60.0
          TROP_Avg_life_of_building_yr = 20.0
          TROP_Biomass_locked_in_construction_material_in_1850_GtBiomass = 30.0
          TROP_clear_cut_fraction = 0.5
          TROP_Dead_biomass__litter_and_soil_organic_matter_SOM_in_1850_GtBiomass = 160.0
          TROP_DeadB_and_SOM_densitiy_in_1850_tBiomass_pr_km2 = 8500.0
          TROP_Fraction_of_construction_waste_burned_0_1 = 0.5
          TROP_fraction_of_DeadB_and_SOM_being_destroyed_by_clear_cutting = 0.5
          TROP_fraction_of_DeadB_and_SOM_being_destroyed_by_deforesting = 1.0
          TROP_fraction_of_DeadB_and_SOM_being_destroyed_by_energy_harvesting = 0.1
          TROP_fraction_of_DeadB_and_SOM_destroyed_by_natural_fires = 0.0
          TROP_living_biomass_densitiy_in_1850_tBiomass_pr_km2 = 16500.0
          TROP_Living_biomass_in_1850_GtBiomass = 370.0
          TROP_Normal_fire_incidence_fraction_yr = 0.3
          TROP_Ref_historical_deforestation_pct_yr = 1.0
          TROP_Slope_temp_eff_on_potential_biomass_per_km2 = -0.5
          TROP_Speed_of_regrowth_yr = 3.0
          TUNDRA_area_burned_in_1850_Mkm2 = 2.0
          TUNDRA_area_deforested_in_1850_Mkm2 = 0.0
          TUNDRA_area_harvested_in_1850_Mkm2 = 2.5
          TUNDRA_Avg_life_biomass_yr = 100.0
          TUNDRA_Avg_life_of_building_yr = 10.0
          TUNDRA_Biomass_locked_in_construction_material_in_1850_GtBiomass = 1.5
          TUNDRA_Dead_biomass__litter_and_soil_organic_matter_SOM_in_1850_GtBiomass = 1200.0
          TUNDRA_DeadB_and_SOM_densitiy_in_1850_tBiomass_pr_km2 = 65000.0
          TUNDRA_Fraction_of_construction_waste_burned_0_1 = 0.5
          TUNDRA_fraction_of_DeadB_and_SOM_being_destroyed_by_deforesting = 1.0
          TUNDRA_fraction_of_DeadB_and_SOM_being_destroyed_by_energy_harvesting = 0.1
          TUNDRA_fraction_of_DeadB_and_SOM_destroyed_by_natural_fires = 0.0
          TUNDRA_living_biomass_densitiy_in_1850_tBiomass_pr_km2 = 14500.0
          TUNDRA_Living_biomass_in_1850_GtBiomass = 300.0
          TUNDRA_Normal_fire_incidence_fraction_yr = 1.0
          TUNDRA_Ref_historical_deforestation_pct_yr = 0.0
          TUNDRA_runoff_time = 2000.0
          TUNDRA_Speed_of_regrowth_yr = 3.0
          TUNDRA_Time_to_decompose_undisturbed_dead_biomass_yr = 1000.0
          UNIT_conversion_1_km3 = 1.0
          UNIT_conversion_1_yr = 1.0
          UNIT_conversion_C_to_pH = 1.0
          UNIT_Conversion_from__km3__km_yr___to_Mkm2_yr = 1.0e-6
          UNIT_conversion_from_km_to_m = 1000.0
          UNIT_Conversion_from_km3_to_km2 = 1.0
          UNIT_Conversion_from_N2O_amount_to_concentration_ppb_MtN2O = 0.305
          UNIT_conversion_Gm3_to_km3 = 1.0
          UNIT_conversion_Gt_to_kt = 1.0e6
          UNIT_conversion_Gt_to_Mt = 1000.0
          UNIT_conversion_GtBiomass_yr_to_Mkm2_yr = 1000.0
          UNIT_conversion_GtC_to_MtC = 1000.0
          UNIT_conversion_GtIce_to_ZJ_melting = 1.0
          UNIT_conversion_km2___km_to_km3 = 1.0
          UNIT_conversion_km2_to_Mkm2 = 1.0e6
          UNIT_conversion_km3_to_Gm3 = 1.0
          UNIT_conversion_km3_km_to_km2 = 1.0
          UNIT_conversion_m2_to_km2 = 1.0e6
          UNIT_conversion_m2_to_Mkm2 = 1.0e12
          UNIT_conversion_Sv_to_Gm3_yr = 31536.0
          UNIT_conversion_to_Gm3 = 1.0
          UNIT_conversion_to_km2_yr = 1.0
          UNIT_conversion_to_yr = 1.0
          UNIT_conversion_W_to_ZJ_s = 1.0
          UNIT_conversion_ymoles___litre_to_dless = 1.0
          UNIT_conversion_yr_to_dless = 1.0
          Urban_area_fraction_2000 = 0.004
          Use_of_GRASS_biomass_for_construction_in_1850_pct = 0.05
          Use_of_GRASS_biomass_for_energy_in_1850_pct = 1.0
          Use_of_NF_biomass_for_construction_in_1850_pct = 0.58
          Use_of_NF_biomass_for_energy_in_1850_pct = 1.09
          Use_of_TROP_biomass_for_construction_in_1850_pct = 0.48
          Use_of_TROP_biomass_for_energy_in_1850_pct = 0.07
          Use_of_TUNDRA_biomass_for_construction_in_1850_pct = 0.05
          Use_of_TUNDRA_biomass_for_energy_in_1850_pct = 1.0
          VAES_puls_repetition = 40.0
          VAES_pulse_duration = 10.0
          VAES_pulse_height = 1.0
          Value_of_anthropogenic_aerosol_emissions_during_2015 = 0.225
          Water_content_of_evaporation_g_kg_per_ZJ_yr = 0.00125
          Wetlands_area_1850 = 1.0e7
          When_first_destroyed_yr = float(2020)
          When_methane_hydrates_first_released_yr = float(2020)
          When_to_sample_for_CO2_experiment_yr = float(20000000)
          Yr_to_cut_mm_emi_abrubtly_to_zero_y = 2020.0
          Zero_C_on_K_scale_K = 273.15
          Zetta = 1.0e21
          CO2_concentration_in_1750_ppm = 2.0
          N2O_ie_N_1750_ppb = 2.0
          CH4_ie_M_1750_ppb = 2.0
          LW_Clear_sky_emissions_from_atm_W_m2_in_1850 = 2.0
          SW_surface_absorption_W_m2_in_1850 = 2.0
          SW_surface_reflection_W_m2_in_1850 = 2.0
          C_in_TUNDRA_DeadB_and_soil_in_1850_GtC = 2.0
          C_in_TUNDRA_LB_in_1850_GtC = 2.0
          Ga__BB_radiation_less_TOA_radiation_W_m2_in_1850 = 2.0
          Biomass_new_growing_1850_GtBiomass___yr = 2.0
          Switch_to_choose_CO2_CH4_aerosol_other_GHG_0normal_1zero_from_202constant_from_2010 = 2.0
          LW_TOA_radiation_from_atm_to_space_in_1850_W_m2 = 2.0
          equationConstructors = Function[]
        end
        function generateEquations0()
          println("#Equation generated:" * "50" * "in: " * "generateEquations0")
          [
            0 ~ var"combi_E3_SC_1_CO2_GtC_yr_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_E3_SC_1_CO2_GtC_yr_tableID, 1, combi_E3_SC_1_CO2_GtC_yr_u),
            0 ~ var"combi_E3_SC_1_CH4_GtC_yr_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_E3_SC_1_CH4_GtC_yr_tableID, 1, combi_E3_SC_1_CH4_GtC_yr_u),
            0 ~ var"combi_E3_SC_1_N2O_Mt_yr_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_E3_SC_1_N2O_Mt_yr_tableID, 1, combi_E3_SC_1_N2O_Mt_yr_u),
            0 ~
              var"combi_E3_SC_1_Kyoto_F_kt_yr_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_E3_SC_1_Kyoto_F_kt_yr_tableID, 1, combi_E3_SC_1_Kyoto_F_kt_yr_u),
            0 ~
              var"combi_E3_SC_1_Montreal_gases_kt_yr_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_E3_SC_1_Montreal_gases_kt_yr_tableID, 1, combi_E3_SC_1_Montreal_gases_kt_yr_u),
            0 ~ var"combi_E3_SC_2_CO2_GtC_yr_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_E3_SC_2_CO2_GtC_yr_tableID, 1, combi_E3_SC_2_CO2_GtC_yr_u),
            0 ~ var"combi_E3_SC_2_CH4_GtC_yr_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_E3_SC_2_CH4_GtC_yr_tableID, 1, combi_E3_SC_2_CH4_GtC_yr_u),
            0 ~ var"combi_E3_SC_2_N2O_Mt_yr_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_E3_SC_2_N2O_Mt_yr_tableID, 1, combi_E3_SC_2_N2O_Mt_yr_u),
            0 ~
              var"combi_E3_SC_2_Kyoto_F_kt_yr_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_E3_SC_2_Kyoto_F_kt_yr_tableID, 1, combi_E3_SC_2_Kyoto_F_kt_yr_u),
            0 ~
              var"combi_E3_SC_2_Montreal_gases_kt_yr_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_E3_SC_2_Montreal_gases_kt_yr_tableID, 1, combi_E3_SC_2_Montreal_gases_kt_yr_u),
            0 ~ var"combi_E3_SC_3_CO2_GtC_yr_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_E3_SC_3_CO2_GtC_yr_tableID, 1, combi_E3_SC_3_CO2_GtC_yr_u),
            0 ~ var"combi_E3_SC_3_CH4_GtC_yr_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_E3_SC_3_CH4_GtC_yr_tableID, 1, combi_E3_SC_3_CH4_GtC_yr_u),
            0 ~ var"combi_E3_SC_3_N2O_Mt_yr_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_E3_SC_3_N2O_Mt_yr_tableID, 1, combi_E3_SC_3_N2O_Mt_yr_u),
            0 ~
              var"combi_E3_SC_3_Kyoto_F_kt_yr_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_E3_SC_3_Kyoto_F_kt_yr_tableID, 1, combi_E3_SC_3_Kyoto_F_kt_yr_u),
            0 ~
              var"combi_E3_SC_3_Montreal_gases_kt_yr_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_E3_SC_3_Montreal_gases_kt_yr_tableID, 1, combi_E3_SC_3_Montreal_gases_kt_yr_u),
            0 ~ var"combi_E3_SC_4_CO2_GtC_yr_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_E3_SC_4_CO2_GtC_yr_tableID, 1, combi_E3_SC_4_CO2_GtC_yr_u),
            0 ~ var"combi_E3_SC_4_CH4_GtC_yr_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_E3_SC_4_CH4_GtC_yr_tableID, 1, combi_E3_SC_4_CH4_GtC_yr_u),
            0 ~ var"combi_E3_SC_4_N2O_Mt_yr_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_E3_SC_4_N2O_Mt_yr_tableID, 1, combi_E3_SC_4_N2O_Mt_yr_u),
            0 ~
              var"combi_E3_SC_4_Kyoto_F_kt_yr_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_E3_SC_4_Kyoto_F_kt_yr_tableID, 1, combi_E3_SC_4_Kyoto_F_kt_yr_u),
            0 ~
              var"combi_E3_SC_4_Montreal_gases_kt_yr_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_E3_SC_4_Montreal_gases_kt_yr_tableID, 1, combi_E3_SC_4_Montreal_gases_kt_yr_u),
            0 ~
              var"combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_MtCH4_yr_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_MtCH4_yr_tableID,
                1,
                combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_MtCH4_yr_u,
              ),
            0 ~
              var"combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP3_MtCH4_yr_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP3_MtCH4_yr_tableID,
                1,
                combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP3_MtCH4_yr_u,
              ),
            0 ~
              var"combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP45_MtCH4_yr_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP45_MtCH4_yr_tableID,
                1,
                combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP45_MtCH4_yr_u,
              ),
            0 ~
              var"combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP6_MtCH4_yr_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP6_MtCH4_yr_tableID,
                1,
                combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP6_MtCH4_yr_u,
              ),
            0 ~
              var"combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP85_MtCH4_yr_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP85_MtCH4_yr_tableID,
                1,
                combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP85_MtCH4_yr_u,
              ),
            0 ~
              var"combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_GtC_yr_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_GtC_yr_tableID,
                1,
                combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_GtC_yr_u,
              ),
            0 ~
              var"combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_GtC_yr_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_GtC_yr_tableID,
                1,
                combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_GtC_yr_u,
              ),
            0 ~
              var"combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_hist_with_JR__large_scale_CCS__exp_GtC_yr_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_hist_with_JR__large_scale_CCS__exp_GtC_yr_tableID,
                1,
                combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_hist_with_JR__large_scale_CCS__exp_GtC_yr_u,
              ),
            0 ~
              var"combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP3_GtC_yr_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP3_GtC_yr_tableID,
                1,
                combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP3_GtC_yr_u,
              ),
            0 ~
              var"combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP45_GtC_yr_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP45_GtC_yr_tableID,
                1,
                combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP45_GtC_yr_u,
              ),
            0 ~
              var"combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP6_GtC_yr_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP6_GtC_yr_tableID,
                1,
                combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP6_GtC_yr_u,
              ),
            0 ~
              var"combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP85_GtC_yr_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP85_GtC_yr_tableID,
                1,
                combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP85_GtC_yr_u,
              ),
            0 ~
              var"combi_Emissions_of_anthro_N20_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_MtN2O_yr_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_Emissions_of_anthro_N20_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_MtN2O_yr_tableID,
                1,
                combi_Emissions_of_anthro_N20_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_MtN2O_yr_u,
              ),
            0 ~
              var"combi_Emissions_of_anthro_N20_with_JR_2052_shape_MtN20_yr_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_Emissions_of_anthro_N20_with_JR_2052_shape_MtN20_yr_tableID,
                1,
                combi_Emissions_of_anthro_N20_with_JR_2052_shape_MtN20_yr_u,
              ),
            0 ~
              var"combi_Emissions_of_Kyoto_Flour_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_kt_yr_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_Emissions_of_Kyoto_Flour_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_kt_yr_tableID,
                1,
                combi_Emissions_of_Kyoto_Flour_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_kt_yr_u,
              ),
            0 ~
              var"combi_Emissions_of_Kyoto_Flour_with_JR_2052_shape_kt_yr_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_Emissions_of_Kyoto_Flour_with_JR_2052_shape_kt_yr_tableID,
                1,
                combi_Emissions_of_Kyoto_Flour_with_JR_2052_shape_kt_yr_u,
              ),
            0 ~
              var"combi_Emissions_of_Montreal_gases_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_kt_yr_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_Emissions_of_Montreal_gases_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_kt_yr_tableID,
                1,
                combi_Emissions_of_Montreal_gases_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_kt_yr_u,
              ),
            0 ~
              var"combi_Emissions_of_Montreal_gases_with_JR_2052_shape_kt_yr_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_Emissions_of_Montreal_gases_with_JR_2052_shape_kt_yr_tableID,
                1,
                combi_Emissions_of_Montreal_gases_with_JR_2052_shape_kt_yr_u,
              ),
            0 ~
              var"combi_CH4_emissions_from_CO2e_C_Roads_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_CH4_emissions_from_CO2e_C_Roads_tableID, 1, combi_CH4_emissions_from_CO2e_C_Roads_u),
            0 ~
              var"combi_CH4_emissions_from_CO2e_CAT_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_CH4_emissions_from_CO2e_CAT_tableID, 1, combi_CH4_emissions_from_CO2e_CAT_u),
            0 ~
              var"combi_CH4_emissions_pct_contribution_to_Total_CO2e_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_CH4_emissions_pct_contribution_to_Total_CO2e_tableID,
                1,
                combi_CH4_emissions_pct_contribution_to_Total_CO2e_u,
              ),
            0 ~
              var"combi_CO2_emissions_from_CO2e_C_Roads_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_CO2_emissions_from_CO2e_C_Roads_tableID, 1, combi_CO2_emissions_from_CO2e_C_Roads_u),
            0 ~
              var"combi_CO2_emissions_from_CO2e_CAT_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_CO2_emissions_from_CO2e_CAT_tableID, 1, combi_CO2_emissions_from_CO2e_CAT_u),
            0 ~
              var"combi_CO2_emissions_pct_contribution_to_Total_CO2e_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_CO2_emissions_pct_contribution_to_Total_CO2e_tableID,
                1,
                combi_CO2_emissions_pct_contribution_to_Total_CO2e_u,
              ),
            0 ~
              var"combi_Historical_aerosol_emissions_anthro_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_Historical_aerosol_emissions_anthro_tableID, 1, combi_Historical_aerosol_emissions_anthro_u),
            0 ~
              var"combi_Historical_forcing_from_solar_insolation_W_m2_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_Historical_forcing_from_solar_insolation_W_m2_tableID,
                1,
                combi_Historical_forcing_from_solar_insolation_W_m2_u,
              ),
            0 ~
              var"combi_Historical_aerosol_forcing_volcanic_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_Historical_aerosol_forcing_volcanic_tableID, 1, combi_Historical_aerosol_forcing_volcanic_u),
            0 ~
              var"combi_OGHG_Kyoto_Flour_emi_rcp3_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_OGHG_Kyoto_Flour_emi_rcp3_tableID, 1, combi_OGHG_Kyoto_Flour_emi_rcp3_u),
            0 ~
              var"combi_OGHG_Kyoto_Flour_emi_rcp45_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_OGHG_Kyoto_Flour_emi_rcp45_tableID, 1, combi_OGHG_Kyoto_Flour_emi_rcp45_u),
            0 ~
              var"combi_OGHG_Kyoto_Flour_emi_rcp6_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_OGHG_Kyoto_Flour_emi_rcp6_tableID, 1, combi_OGHG_Kyoto_Flour_emi_rcp6_u),
          ]
        end
        function generateEquations1()
          println("#Equation generated:" * "50" * "in: " * "generateEquations1")
          [
            0 ~
              var"combi_OGHG_Kyoto_Flour_emi_rcp85_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_OGHG_Kyoto_Flour_emi_rcp85_tableID, 1, combi_OGHG_Kyoto_Flour_emi_rcp85_u),
            0 ~
              var"combi_Kyoto_Flour_emissions_from_CO2e_C_Roads_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_Kyoto_Flour_emissions_from_CO2e_C_Roads_tableID, 1, combi_Kyoto_Flour_emissions_from_CO2e_C_Roads_u),
            0 ~
              var"combi_Kyoto_Flour_emissions_from_CO2e_CAT_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_Kyoto_Flour_emissions_from_CO2e_CAT_tableID, 1, combi_Kyoto_Flour_emissions_from_CO2e_CAT_u),
            0 ~
              var"combi_Kyoto_Flour_emissions_pct_contribution_to_Total_CO2e_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_Kyoto_Flour_emissions_pct_contribution_to_Total_CO2e_tableID,
                1,
                combi_Kyoto_Flour_emissions_pct_contribution_to_Total_CO2e_u,
              ),
            0 ~
              var"combi_OGHG_Montreal_gases_emi_rcp3_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_OGHG_Montreal_gases_emi_rcp3_tableID, 1, combi_OGHG_Montreal_gases_emi_rcp3_u),
            0 ~
              var"combi_OGHG_Montreal_gases_emi_rcp45_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_OGHG_Montreal_gases_emi_rcp45_tableID, 1, combi_OGHG_Montreal_gases_emi_rcp45_u),
            0 ~
              var"combi_OGHG_Montreal_gases_emi_rcp6_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_OGHG_Montreal_gases_emi_rcp6_tableID, 1, combi_OGHG_Montreal_gases_emi_rcp6_u),
            0 ~
              var"combi_OGHG_Montreal_gases_emi_rcp85_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_OGHG_Montreal_gases_emi_rcp85_tableID, 1, combi_OGHG_Montreal_gases_emi_rcp85_u),
            0 ~
              var"combi_othGHG_N20_man_made_emissions_rcp3_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_othGHG_N20_man_made_emissions_rcp3_tableID, 1, combi_othGHG_N20_man_made_emissions_rcp3_u),
            0 ~
              var"combi_othGHG_N20_man_made_emissions_rcp45_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_othGHG_N20_man_made_emissions_rcp45_tableID, 1, combi_othGHG_N20_man_made_emissions_rcp45_u),
            0 ~
              var"combi_othGHG_N20_man_made_emissions_rcp6_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_othGHG_N20_man_made_emissions_rcp6_tableID, 1, combi_othGHG_N20_man_made_emissions_rcp6_u),
            0 ~
              var"combi_othGHG_N20_man_made_emissions_rcp85_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_othGHG_N20_man_made_emissions_rcp85_tableID, 1, combi_othGHG_N20_man_made_emissions_rcp85_u),
            0 ~
              var"combi_RCP_3_CO2_concentration_1850_2100_ppm_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_RCP_3_CO2_concentration_1850_2100_ppm_tableID, 1, combi_RCP_3_CO2_concentration_1850_2100_ppm_u),
            0 ~
              var"combi_RCP_45_CO2_concentration_1850_2100_ppm_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_RCP_45_CO2_concentration_1850_2100_ppm_tableID, 1, combi_RCP_45_CO2_concentration_1850_2100_ppm_u),
            0 ~
              var"combi_RCP_6_CO2_concentration_1850_2100_ppm_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_RCP_6_CO2_concentration_1850_2100_ppm_tableID, 1, combi_RCP_6_CO2_concentration_1850_2100_ppm_u),
            0 ~
              var"combi_RCP_85_CO2_concentration_1850_2100_ppm_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_RCP_85_CO2_concentration_1850_2100_ppm_tableID, 1, combi_RCP_85_CO2_concentration_1850_2100_ppm_u),
            0 ~
              var"combi_Montreal_gases_emissions_from_CO2e_C_Roads_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_Montreal_gases_emissions_from_CO2e_C_Roads_tableID,
                1,
                combi_Montreal_gases_emissions_from_CO2e_C_Roads_u,
              ),
            0 ~
              var"combi_Montreal_gases_emissions_from_CO2e_CAT_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_Montreal_gases_emissions_from_CO2e_CAT_tableID, 1, combi_Montreal_gases_emissions_from_CO2e_CAT_u),
            0 ~
              var"combi_Montreal_gases_emissions_pct_contribution_to_Total_CO2e_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_Montreal_gases_emissions_pct_contribution_to_Total_CO2e_tableID,
                1,
                combi_Montreal_gases_emissions_pct_contribution_to_Total_CO2e_u,
              ),
            0 ~
              var"combi_N2O_man_made_emissions_from_CO2e_C_Roads_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_N2O_man_made_emissions_from_CO2e_C_Roads_tableID,
                1,
                combi_N2O_man_made_emissions_from_CO2e_C_Roads_u,
              ),
            0 ~
              var"combi_N2O_man_made_emissions_from_CO2e_CAT_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_N2O_man_made_emissions_from_CO2e_CAT_tableID, 1, combi_N2O_man_made_emissions_from_CO2e_CAT_u),
            0 ~
              var"combi_N2O_emissions_pct_contribution_to_Total_CO2e_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_N2O_emissions_pct_contribution_to_Total_CO2e_tableID,
                1,
                combi_N2O_emissions_pct_contribution_to_Total_CO2e_u,
              ),
            0 ~
              var"combi_Sea_level_rise_history_mm_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_Sea_level_rise_history_mm_tableID, 1, combi_Sea_level_rise_history_mm_u),
            0 ~
              var"combi_All_Human_activity_emissions_GtCO2e_yr_Base_for_tipping_point_search_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_All_Human_activity_emissions_GtCO2e_yr_Base_for_tipping_point_search_tableID,
                1,
                combi_All_Human_activity_emissions_GtCO2e_yr_Base_for_tipping_point_search_u,
              ),
            0 ~
              var"combi_Arctic_freezing_cutoff_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_Arctic_freezing_cutoff_tableID, 1, combi_Arctic_freezing_cutoff_u),
            0 ~
              var"combi_Blocked_by_H20_hist_Table_lookup_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_Blocked_by_H20_hist_Table_lookup_tableID, 1, combi_Blocked_by_H20_hist_Table_lookup_u),
            0 ~
              var"combi_Blocked_by_H20_Table_lookup_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_Blocked_by_H20_Table_lookup_tableID, 1, combi_Blocked_by_H20_Table_lookup_u),
            0 ~
              var"combi_Effect_of_heat_in_atm_on_melting_ice__cut_off__y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_Effect_of_heat_in_atm_on_melting_ice__cut_off__tableID,
                1,
                combi_Effect_of_heat_in_atm_on_melting_ice__cut_off__u,
              ),
            0 ~
              var"combi_Exp_12a_reduction_in_emissions_LOOKUP_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_Exp_12a_reduction_in_emissions_LOOKUP_tableID, 1, combi_Exp_12a_reduction_in_emissions_LOOKUP_u),
            0 ~
              var"combi_EXP_12b_CCS_from_2015_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_EXP_12b_CCS_from_2015_tableID, 1, combi_EXP_12b_CCS_from_2015_u),
            0 ~
              var"combi_EXP_12e_white_surfaces_ease_in_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_EXP_12e_white_surfaces_ease_in_tableID, 1, combi_EXP_12e_white_surfaces_ease_in_u),
            0 ~
              var"combi_Fraction_blocked_by_CH4_spectrum_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_Fraction_blocked_by_CH4_spectrum_tableID, 1, combi_Fraction_blocked_by_CH4_spectrum_u),
            0 ~
              var"combi_Fraction_blocked_by_CO2_spectrum_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_Fraction_blocked_by_CO2_spectrum_tableID, 1, combi_Fraction_blocked_by_CO2_spectrum_u),
            0 ~
              var"combi_Future_shape_of_anthropogenic_aerosol_emissions_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_Future_shape_of_anthropogenic_aerosol_emissions_tableID,
                1,
                combi_Future_shape_of_anthropogenic_aerosol_emissions_u,
              ),
            0 ~
              var"combi_Melting_constraint_from_the_heat_in__ocean__surface_reservoir_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_Melting_constraint_from_the_heat_in__ocean__surface_reservoir_tableID,
                1,
                combi_Melting_constraint_from_the_heat_in__ocean__surface_reservoir_u,
              ),
            0 ~
              var"combi_Melting_constraint_from_the_heat_in_atmosphere_reservoir_fraction_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_Melting_constraint_from_the_heat_in_atmosphere_reservoir_fraction_tableID,
                1,
                combi_Melting_constraint_from_the_heat_in_atmosphere_reservoir_fraction_u,
              ),
            0 ~
              var"combi_Melting_restraint_for_permafrost_from_heat_in_atmophere_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_Melting_restraint_for_permafrost_from_heat_in_atmophere_tableID,
                1,
                combi_Melting_restraint_for_permafrost_from_heat_in_atmophere_u,
              ),
            0 ~
              var"combi_NATURE_CCS_removal_experiment_multiplier_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(
                combi_NATURE_CCS_removal_experiment_multiplier_tableID,
                1,
                combi_NATURE_CCS_removal_experiment_multiplier_u,
              ),
            0 ~
              var"combi_NF_clear_cut_fraction_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_NF_clear_cut_fraction_tableID, 1, combi_NF_clear_cut_fraction_u),
            0 ~ var"combi_NF_usage_cutoff_y[1]" - Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_NF_usage_cutoff_tableID, 1, combi_NF_usage_cutoff_u),
            0 ~
              var"combi_Permafrost_melting_cutoff_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_Permafrost_melting_cutoff_tableID, 1, combi_Permafrost_melting_cutoff_u),
            0 ~
              var"combi_RCPFossil_fuel_usage_cutoff_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_RCPFossil_fuel_usage_cutoff_tableID, 1, combi_RCPFossil_fuel_usage_cutoff_u),
            0 ~
              var"combi_Snowball_earth_cutoff_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_Snowball_earth_cutoff_tableID, 1, combi_Snowball_earth_cutoff_u),
            0 ~
              var"combi_Thermal_expansion_deep_in_1850_pct_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_Thermal_expansion_deep_in_1850_pct_tableID, 1, combi_Thermal_expansion_deep_in_1850_pct_u),
            0 ~
              var"combi_Thermal_expansion_deep_pct_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_Thermal_expansion_deep_pct_tableID, 1, combi_Thermal_expansion_deep_pct_u),
            0 ~
              var"combi_Thermal_expansion_surface_in_1850_pct_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_Thermal_expansion_surface_in_1850_pct_tableID, 1, combi_Thermal_expansion_surface_in_1850_pct_u),
            0 ~
              var"combi_Thermal_expansion_surface_pct_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_Thermal_expansion_surface_pct_tableID, 1, combi_Thermal_expansion_surface_pct_u),
            0 ~
              var"combi_TROP_deforestation_cutoff_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_TROP_deforestation_cutoff_tableID, 1, combi_TROP_deforestation_cutoff_u),
            0 ~
              var"combi_TROP_deforestation_cutoff_effect_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_TROP_deforestation_cutoff_effect_tableID, 1, combi_TROP_deforestation_cutoff_effect_u),
            0 ~
              var"combi_TROP_deforestion_multiplier_wrt_2000_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_TROP_deforestion_multiplier_wrt_2000_tableID, 1, combi_TROP_deforestion_multiplier_wrt_2000_u),
          ]
        end
        function generateEquations2()
          println("#Equation generated:" * "50" * "in: " * "generateEquations2")
          [
            0 ~
              var"combi_Urbanzation_Effect_on_biomass_use_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_Urbanzation_Effect_on_biomass_use_tableID, 1, combi_Urbanzation_Effect_on_biomass_use_u),
            0 ~
              var"combi_Population_Lookup_bn_y[1]" -
              Modelica_Blocks_Tables_Internal_getTable1DValueNoDer2(combi_Population_Lookup_bn_tableID, 1, combi_Population_Lookup_bn_u),
            0 ~ Time - t,
            D(Antarctic_ice_volume_km3) ~ -flow_Antarctic_ice_melting__pos__or_freezing__neg__km3_yr,
            D(Arctic_ice__on_sea__area_km2) ~ -flow_Arctic_ice_melting__pos__or_freezing__neg__km2_yr,
            D(C_in_atmosphere_GtC) ~
              (
                (
                  (
                    (
                      (
                        (
                          (
                            (
                              flow_Avg_volcanic_activity_GtC_yr +
                              flow_CH4_in_the_atmosphere_converted_to_CO2 +
                              flow_CO2_flux_GRASS_to_atm_Gtc_yr +
                              flow_CO2_flux_NF_to_atm_Gtc_yr +
                              flow_CO2_flux_TROP_to_atm_GtC_yr +
                              flow_CO2_flux_TUNDRA_to_atm_Gtc_yr +
                              flow_C_release_from_permafrost_melting_as_CO2_GtC_yr +
                              flow_Man_made_fossil_C_emissions_GtC_yr +
                              flow_Methane_hydrates_released_and_converted_to_CO2_by_bacteria_GtC
                            ) - flow_CO2_flux_from_atm_to_GRASS_for_new_growth_GtC_yr
                          ) - flow_CO2_flux_from_atm_to_NF_for_new_growth_GtC_yr
                        ) - flow_CO2_flux_from_atm_to_TROP_for_new_growth_GtC_yr
                      ) - flow_CO2_flux_from_atm_to_TUNDRA_for_new_growth_GtC_yr
                    ) - flow_C_diffusion_into_ocean_from_atm
                  ) - flow_C_removal_rate_from_atm_for_nature_May_2020_GtC_y
                ) - flow_Carbon_captured_and_stored_GtC___yr
              ) - flow_NATURE_CCS_Fig3_GtC_yr,
            D(C_in_atmosphere_in_form_of_CH4) ~
              (flow_CH4_release___capture_from_permafrost_area_loss___gain_GtC_yr + flow_Human_activity_CH4_emissions + flow_Methanehydrate_experimental_release_GtC__yr + flow_Natural_CH4_emissions) -
              flow_CH4_conversion_to_CO2_and_H2O,
            D(C_in_cold_surface_water_GtC) ~ (flow_C_diffusion_into_ocean_from_atm + flow_Carbon_flow_from_warm_to_cold_surface_GtC_per_yr) - flow_Carbon_flow_from_cold_surface_downwelling_Gtc_per_yr,
            D(C_in_cold_water_trunk_downwelling_GtC) ~ flow_Carbon_flow_from_cold_surface_downwelling_Gtc_per_yr - flow_Carbon_flow_from_cold_to_deep_GtC_per_yr,
            D(C_in_deep_water_volume_1km_to_bottom_GtC) ~ (flow_Carbon_flow_from_cold_to_deep_GtC_per_yr - flow_Carbon_flow_from_deep) - flow_Depositing_of_C_to_sediment,
            D(C_in_intermediate_upwelling_water_100m_to_1km_GtC) ~ flow_Carbon_flow_from_deep - flow_Carbon_flow_from_intermediate_to_surface_box_GtC_per_yr,
            D(C_in_permafrost_in_form_of_CH4) ~ -flow_CH4_release___capture_from_permafrost_area_loss___gain_GtC_yr,
            D(C_in_sediment) ~ flow_Biological_removal_of_C_from_WSW_GtC_per_yr + flow_Depositing_of_C_to_sediment,
            D(C_in_warm_surface_water_GtC) ~
              ((flow_C_runoff_from_biomass_soil + flow_Carbon_flow_from_intermediate_to_surface_box_GtC_per_yr) - flow_Biological_removal_of_C_from_WSW_GtC_per_yr) -
              flow_Carbon_flow_from_warm_to_cold_surface_GtC_per_yr,
            D(Cold_surface_water_volume_Gm3) ~ flow_Flow_of_water_from_warm_to_cold_surface_Gcubicm_per_yr - flow_Flow_of_cold_surface_water_welling_down_GcubicM_per_yr,
            D(Cold_water_volume_downwelling_Gm3) ~ flow_Flow_of_cold_surface_water_welling_down_GcubicM_per_yr - flow_Flow_of_cold_water_sinking_to_very_bottom_GcubicM_per_yr,
            D(Cumulative_antarctic_ice_volume_loss_GtIce) ~ flow_Antarctic_ice_losing__pos__or_gaining__neg__GtIce_yr,
            D(Cumulative_C_released_from_permafrost_as_either_CH4_or_CO2_in_GtC_yr) ~ flow_C_released_from_permafrost_as_either_CH4_or_CO2_in_GtC_yr,
            D(Cumulative_carbon_captured_and_stored_GtC) ~ flow_Carbon_captured_and_stored_GtC___yr,
            D(Cumulative_carbon_removed_from_atm_for_nature_May_2020) ~ flow_C_removal_rate_from_atm_for_nature_May_2020_GtC_y,
            D(Cumulative_flow_of_C_to_biomass_since_1850_GtC) ~ flow_Annual_flux_of_C_to_biomass_GtC_pr_yr,
            D(Cumulative_glacial_ice_volume_loss_GtIce) ~ flow_Annual_glacial_ice_losing__pos__or_gaining__neg__GtIce_yr,
            D(Cumulative_Greenland_ice_volume_loss_GtIce) ~ flow_Greenland_ice_losing__pos__or_gaining__neg__GtIce_yr,
            D(Cumulative_heat_to_atm_ZJ) ~ flow_Flow_of_heat_to_atm_ZJ_yr,
            D(Cumulative_ocean_volume_increase_due_to_ice_melting_km3) ~
              flow_Antarctic_ice_melting_as_water_km3_yr +
              flow_Glacial_ice_melting_as_water_km3_yr +
              flow_Greenland_ice_melting_as_water_km3_yr +
              flow_Greenland_ice_melting_that_slid_into_the_ocean_km3_yr,
            D(Cumulative_release_of_C_from_permafrost_GtC) ~ flow_Annual_release_of_C_from_permafrost_GtC_y,
            D(Deep_water_volume_1km_to_4km_Gm3) ~ flow_Flow_of_cold_water_sinking_to_very_bottom_GcubicM_per_yr - flow_Upwelling_from_deep,
            D(DESERT_Mkm2) ~ flow_Shifting_GRASS_to_DESERT_Mkm2_yr - flow_Sifting_DESERT_to_GRASS_Mkm2_yr,
            D(Fossil_fuel_reserves_in_ground_GtC) ~ -flow_Man_made_fossil_C_emissions_GtC_yr,
            D(Glacial_ice_volume_km3) ~ -flow_Glacial_ice_melting__pos__or_freezing__neg__km3_yr,
            D(GRASS_area_burnt_Mkm2) ~ flow_GRASS_burning_Mkm2_yr - flow_GRASS_regrowing_after_being_burnt_Mkm2_yr,
            D(GRASS_area_harvested_Mkm2) ~ flow_GRASS_being_harvested_Mkm2_yr - flow_GRASS_regrowing_after_harvesting_Mkm2_yr,
            D(GRASS_Biomass_locked_in_construction_material_GtBiomass) ~
              (flow_GRASS_for_construction_use_GtBiomass_yr - flow_GRASS_Biomass_in_construction_material_being_burnt_GtBiomass_yr) -
              flow_GRASS_Biomass_in_construction_material_left_to_rot_GtBiomass_yr,
            D(GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass) ~
              (
                (
                  (
                    (
                      (flow_GRASS_Biomass_in_construction_material_left_to_rot_GtBiomass_yr + flow_GRASS_Living_biomass_rotting_GtBiomass_yr) -
                      flow_GRASS_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr
                    ) - flow_GRASS_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr
                  ) - flow_GRASS_Dead_biomass_decomposing_GtBiomass_yr
                ) - flow_GRASS_runoff
              ) - flow_GRASS_soil_degradation_from_forest_fires_GtBiomass_yr,
            D(GRASS_deforested_Mkm2) ~ flow_GRASS_being_deforested_Mkm2_yr - flow_GRASS_regrowing_after_being_deforested_Mkm2_yr,
            D(GRASS_Living_biomass_GtBiomass) ~
              (
                (flow_GRASS_biomass_new_growing_GtBiomass___yr - flow_GRASS_Living_biomass_rotting_GtBiomass_yr) -
                flow_GRASS_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr
              ) - flow_GRASS_for_construction_use_GtBiomass_yr,
            D(GRASS_potential_area_Mkm2) ~
              (
                ((flow_Shifting_NF_to_GRASS_Mkm2_yr + flow_Shifting_TROP_to_GRASS_Mkm2_yr + flow_Sifting_DESERT_to_GRASS_Mkm2_yr) - flow_Shifting_GRASS_to_DESERT_Mkm2_yr) -
                flow_Shifting_GRASS_to_NF_Mkm2_yr
              ) - flow_Shifting_GRASS_to_TROP_Mkm2_yr,
            D(Greenland_ice_volume_on_Greenland_km3) ~ -flow_Greenland_ice_melting__pos__or_freezing__neg__km3_yr - flow_Greenland_ice_sliding_into_the_ocean_km3_yr,
            D(Greenland_ice_volume_that_slid_into_the_ocean_km3) ~
              flow_Greenland_ice_sliding_into_the_ocean_km3_yr - flow_Greenland_ice_that_slid_into_the_ocean_melting__pos__or_freezing__neg__km3_yr,
            D(Heat_in_atmosphere_ZJ) ~
              (
                (
                  (flow_Convection_aka_sensible_heat_flow + flow_Evaporation_aka_latent_heat_flow + flow_LW_surface_emissions_NOT_escaping_through_atm_window + flow_SW_Atmospheric_absorption) -
                  flow_Heat_withdrawn_from_atm_by_melting__pos____added__neg__by_freezing_ice_ZJ_yr
                ) - flow_LW_TOA_radiation_from_atm_to_space
              ) - flow_LW_clear_sky_emissions_to_surface,
            D(Heat_in_deep_ZJ) ~ flow_Heat_flow_from_the_earths_core + flow_Net_heat_flow_ocean_from_surface_to_deep__ZJ_yr_,
            D(Heat_in_surface) ~
              (
                (
                  (
                    (
                      (
                        ((flow_LW_clear_sky_emissions_to_surface + flow_LW_re_radiated_by_clouds + flow_SW_surface_absorption) - flow_Convection_aka_sensible_heat_flow) -
                        flow_Evaporation_aka_latent_heat_flow
                      ) - flow_Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_Greenland_ice_that_slid_into_the_ocean_ZJ_yr
                    ) - flow_Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_antarctic_ice_ZJ_yr
                  ) - flow_Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_arctic_ice_ZJ_yr
                ) - flow_LW_surface_emission
              ) - flow_Net_heat_flow_ocean_from_surface_to_deep__ZJ_yr_,
            D(Intermediate_upwelling_water_volume_100m_to_1km_Gm3) ~ flow_Upwelling_from_deep - flow_Upwelling_to_surface,
            D(Kyoto_Flour_gases_in_atm) ~ flow_Kyoto_Flour_emissions - flow_Kyoto_Flour_degradation,
            D(Montreal_gases_in_atm) ~ flow_Montreal_gases_emissions - flow_Montreal_gases_degradation,
            D(N2O_in_atmosphere_MtN2O) ~ flow_All_N2O_emissions_MtN2O_yr - flow_N2O_degradation_MtN2O_yr,
            D(NATURE_Cumulative_CCS_GtC) ~ flow_NATURE_CCS_Fig3_GtC_yr,
            D(NF_area_burnt_Mkm2) ~ flow_NF_burning_Mkm2_yr - flow_NF_regrowing_after_being_burnt_Mkm2_yr,
            D(NF_area_clear_cut_Mkm2) ~ flow_NF_being_harvested_by_clear_cutting_Mkm2_yr - flow_NF_regrowing_after_being_clear_cut_Mkm2_yr,
            D(NF_area_deforested_Mkm2) ~ flow_NF_being_deforested_Mkm2_yr - flow_NF_regrowing_after_being_deforested_Mkm2_yr,
          ]
        end
        function generateEquations3()
          println("#Equation generated:" * "50" * "in: " * "generateEquations3")
          [
            D(NF_area_harvested_Mkm2) ~ flow_NF_being_harvested_normally_Mkm2_yr - flow_NF_regrowing_after_harvesting_Mkm2_yr,
            D(NF_Biomass_locked_in_construction_material_GtBiomass) ~
              (flow_NF_for_construction_use_GtBiomass_yr - flow_NF_Biomass_in_construction_material_being_burnt_GtBiomass_yr) -
              flow_NF_TUNDRA_Biomass_in_construction_material_left_to_rot_GtBiomass_yr,
            D(NF_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass) ~
              (
                (
                  (
                    (
                      (
                        (flow_NF_Living_biomass_rotting_GtBiomass_yr + flow_NF_TUNDRA_Biomass_in_construction_material_left_to_rot_GtBiomass_yr) -
                        flow_NF_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr
                      ) - flow_NF_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr
                    ) - flow_NF_Dead_biomass_decomposing_GtBiomass_yr
                  ) - flow_NF_runoff
                ) - flow_NF_soil_degradation_from_clear_cutting_GtBiomass_yr
              ) - flow_NF_soil_degradation_from_forest_fires_GtBiomass_yr,
            D(NF_Living_biomass_GtBiomass) ~
              (
                (flow_NF_biomass_new_growing_GtBiomass___yr - flow_NF_Living_biomass_rotting_GtBiomass_yr) -
                flow_NF_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr
              ) - flow_NF_for_construction_use_GtBiomass_yr,
            D(NF_potential_area_Mkm2) ~
              (((flow_Shifting_GRASS_to_NF_Mkm2_yr + flow_Shifting_TROP_to_NF_Mkm2_yr + flow_Shifting_Tundra_to_NF_Mkm2_yr) - flow_Shifting_NF_to_GRASS_Mkm2_yr) - flow_Shifting_NF_to_TROP_Mkm2_yr) -
              flow_Shifting_NF_to_Tundra_Mkm2_yr,
            D(Sum_C_absorbed_by_ocean_GtC) ~ flow_C_absorption_by_ocean_from_atm_for_accumulation,
            D(Sum_heat_to_deep_ocean) ~ flow_Flow_of_heat_to_deep_ocean,
            D(Sum_heat_to_deep_ocean_btw_72_and_08) ~ flow_Flow_of_heat_to_deep_ocean_btw_72_and_08,
            D(Sum_heat_to_surface_ocean_btw_72_and_08) ~ flow_Flow_of_heat_to_surface_ocean_btw_1972_and_2008,
            D(Sum_heat_to_surface_ocean_ZJ) ~ flow_Flow_of_heat_to_surface_ocean,
            D(Sum_man_made_CO2_emissions_GtC) ~ flow_Man_made_fossil_C_emissions_for_cumulation_GtC_yr,
            D(Sum_net_C_to_atm) ~ flow_Net_C_to_atm_rate,
            D(TROP_area_burnt_Mkm2) ~ flow_TROP_burning_Mkm2_yr - flow_TROP_NF_regrowing_after_being_burnt_Mkm2_yr,
            D(TROP_area_clear_cut_Mkm2) ~ flow_TROP_being_harvested_by_clear_cutting_Mkm2_yr - flow_TROP_regrowing_after_being_clear_cut_Mkm2_yr,
            D(TROP_area_deforested_Mkm2) ~ flow_TROP_being_deforested_Mkm2_yr - flow_TROP_regrowing_after_being_deforested_Mkm2_yr,
            D(TROP_area_harvested_Mkm2) ~ flow_TROP_being_harvested_normally_Mkm2_yr - flow_TROP_NF_regrowing_after_harvesting_Mkm2_yr,
            D(TROP_Biomass_locked_in_construction_material_GtBiomass) ~
              (flow_TROP_for_construction_use_GtBiomass_yr - flow_TROP_Biomass_in_construction_material_being_burnt_GtBiomass_yr) - flow_TROP_Biomass_in_construction_material_left_to_rot_GtBiomass_yr,
            D(TROP_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass) ~
              (
                (
                  (
                    (
                      (
                        (flow_TROP_Biomass_in_construction_material_left_to_rot_GtBiomass_yr + flow_TROP_Living_biomass_rotting_GtBiomass_yr) -
                        flow_TROP_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr
                      ) - flow_TROP_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr
                    ) - flow_TROP_Dead_biomass_decomposing_GtBiomass_yr
                  ) - flow_TROP_runoff
                ) - flow_TROP_soil_degradation_from_clear_cutting_GtBiomass_yr
              ) - flow_TROP_soil_degradation_from_forest_fires_GtBiomass_yr,
            D(TROP_Living_biomass_GtBiomass) ~
              (
                (flow_TROP_biomass_new_growing_GtBiomass___yr - flow_TROP_Living_biomass_rotting_GtBiomass_yr) -
                flow_TROP_NF_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr
              ) - flow_TROP_for_construction_use_GtBiomass_yr,
            D(TROP_potential_area_Mkm2) ~ ((flow_Shifting_GRASS_to_TROP_Mkm2_yr + flow_Shifting_NF_to_TROP_Mkm2_yr) - flow_Shifting_TROP_to_GRASS_Mkm2_yr) - flow_Shifting_TROP_to_NF_Mkm2_yr,
            D(TUNDRA_area_burnt_Mkm2) ~ flow_TUNDRA_burning_Mkm2_yr - flow_TUNDRA_regrowing_after_being_burnt_Mkm2_yr,
            D(TUNDRA_area_harvested_Mkm2) ~ flow_TUNDRA_being_harvested_Mkm2_yr - flow_TUNDRA_regrowing_after_harvesting_Mkm2_yr,
            D(TUNDRA_Biomass_locked_in_construction_material_GtBiomass) ~
              (flow_TUNDRA_for_construction_use_GtBiomass_yr - flow_TUNDRA_Biomass_in_construction_material_being_burnt_GtBiomass_yr) -
              flow_TUNDRA_Biomass_in_construction_material_left_to_rot_GtBiomass_yr,
            D(TUNDRA_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass) ~
              (
                (
                  (
                    (
                      (flow_TUNDRA_Biomass_in_construction_material_left_to_rot_GtBiomass_yr + flow_TUNDRA_Living_biomass_rotting_GtBiomass_yr) -
                      flow_TUNDRA_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr
                    ) - flow_TUNDRA_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr
                  ) - flow_TUNDRA_Dead_biomass_decomposing_GtBiomass_yr
                ) - flow_TUNDRA_runoff
              ) - flow_TUNDRA_soil_degradation_from_forest_fires_GtBiomass_yr,
            D(TUNDRA_deforested_Mkm2) ~ flow_TUNDRA_being_deforested_Mkm2_yr - flow_TUNDRA_regrowing_after_being_deforested_Mkm2_yr,
            D(TUNDRA_Living_biomass_GtBiomass) ~
              (
                (flow_TUNDRA_biomass_new_growing_GtBiomass___yr - flow_TUNDRA_Living_biomass_rotting_GtBiomass_yr) -
                flow_TUNDRA_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr
              ) - flow_TUNDRA_for_construction_use_GtBiomass_yr,
            D(Tundra_potential_area_Mkm2) ~
              ((flow_Shifting_NF_to_Tundra_Mkm2_yr + flow_Shifting_ice_on_land_to_tundra_Mkm2_yr) - flow_Shifting_Tundra_to_NF_Mkm2_yr) - flow_Shifting_tundra_to_ice_on_land_Mkm2_yr,
            D(Volcanic_aerosols_in_stratosphere) ~ flow_Volcanic_aerosols_emissions - flow_Volcanic_aerosols_removed_from_stratosphere,
            D(Warm_surface_water_volume_Gm3) ~ flow_Upwelling_to_surface - flow_Flow_of_water_from_warm_to_cold_surface_Gcubicm_per_yr,
            D(Wetlands_area) ~ -flow_Rate_of_destruction_of_wetlands,
            0 ~ combi_All_Human_activity_emissions_GtCO2e_yr_Base_for_tipping_point_search_u - t,
            0 ~ combi_Arctic_freezing_cutoff_u + -Arctic_ice__on_sea__area_km2 / Arctic_ice_area_max_km2,
            0 ~ combi_Blocked_by_H20_hist_Table_lookup_u - Humidity_of_atmosphere_current_g_kg,
            0 ~ combi_Blocked_by_H20_Table_lookup_u - Humidity_of_atmosphere_current_g_kg,
            0 ~ combi_Effect_of_heat_in_atm_on_melting_ice__cut_off__u - 0.0009749724570280889Heat_in_atmosphere_ZJ,
            0 ~ combi_Exp_12a_reduction_in_emissions_LOOKUP_u - t,
            0 ~ combi_EXP_12b_CCS_from_2015_u - t,
            0 ~ combi_EXP_12e_white_surfaces_ease_in_u - t,
            0 ~ combi_Fraction_blocked_by_CH4_spectrum_u - CH4_concentration_ppb,
            0 ~ combi_Fraction_blocked_by_CO2_spectrum_u - CO2_concentration_ppm,
            0 ~ combi_Future_shape_of_anthropogenic_aerosol_emissions_u - t,
            0 ~ combi_Melting_constraint_from_the_heat_in__ocean__surface_reservoir_u - Ocean_heat_used_for_melting_last_year_ZJ_yr,
            0 ~ combi_Melting_constraint_from_the_heat_in_atmosphere_reservoir_fraction_u - Atmos_heat_used_for_melting_last_year_1_yr,
            0 ~ combi_Melting_restraint_for_permafrost_from_heat_in_atmophere_u + -Heat_gained___needed_for_the_desired_freezing___unfreezing_of_permafrost_ZJ_yr / Heat_in_atmosphere_ZJ,
            0 ~ combi_NATURE_CCS_removal_experiment_multiplier_u - t,
            0 ~ combi_NF_clear_cut_fraction_u - t,
            0 ~ combi_NF_usage_cutoff_u - NF_usage_as_pct_of_potial_area,
            0 ~ combi_Permafrost_melting_cutoff_u - C_in_permafrost_in_form_of_CH4,
            0 ~ combi_RCPFossil_fuel_usage_cutoff_u + -Fossil_fuel_reserves_in_ground_GtC / Fossil_fuel_reserves_in_ground_1850_GtC,
            0 ~ combi_Snowball_earth_cutoff_u + -Land_covered_with_ice_km2 / Land_area_km2,
          ]
        end
        function generateEquations4()
          println("#Equation generated:" * "50" * "in: " * "generateEquations4")
          [
            0 ~ combi_Thermal_expansion_deep_in_1850_pct_u - 4.0,
            0 ~ combi_Thermal_expansion_deep_pct_u - Temp__ocean__deep_in_C,
            0 ~ combi_Thermal_expansion_surface_in_1850_pct_u - 13.66500000000002,
            0 ~ combi_Thermal_expansion_surface_pct_u - Temp_surface_C,
            0 ~ combi_TROP_deforestation_cutoff_u - TROP_deforested_as_pct_of_potial_area,
            0 ~ combi_TROP_deforestation_cutoff_effect_u - TROP_deforested_as_pct_of_potial_area,
            0 ~ combi_TROP_deforestion_multiplier_wrt_2000_u - t,
            0 ~ combi_Urbanzation_Effect_on_biomass_use_u - t,
            0 ~ combi_Population_Lookup_bn_u - Time,
            D(Arctic_land_surface_temp_anomaly_compared_to_1850) ~ 0.04 * (Temp_surface_anomaly_compared_to_1850_degC - Arctic_land_surface_temp_anomaly_compared_to_1850),
            D(Biological_removal_of_C_from_WSW_GtC_per_yr) ~ 0.04 * (Net_marine_primary_production_NMPP_GtC_pr_yr - Biological_removal_of_C_from_WSW_GtC_per_yr),
            D(Effect_of_temp_on_permafrost_melting_dmnl) ~
              0.2 * ((1.0 + Slope_btw_temp_and_permafrost_melting___freezing * (0.25Temp_diff_relevant_for_melting_or_freezing_from_1850 - 1.0)) - Effect_of_temp_on_permafrost_melting_dmnl),
            D(Temp_diff_relevant_for_melting_or_freezing_arctic_ice_from_1850) ~
              0.06666666666666667 * (Temp_surface_anomaly_compared_to_1850_degC - Temp_diff_relevant_for_melting_or_freezing_arctic_ice_from_1850),
            D(Temp_diff_relevant_for_melting_or_freezing_from_1850) ~ 0.3333333333333333 * ((Temp_surface_C - 13.66500000000002) - Temp_diff_relevant_for_melting_or_freezing_from_1850),
            D(yr_on_yr_change_in_C_in_atm_GtC_yr) ~ (C_in_atmosphere_GtC - C_in_atm_1_yr_ago_GtC) - yr_on_yr_change_in_C_in_atm_GtC_yr,
            D(C_in_ocean_1_yr_ago_GtC) ~ -((C_in_ocean_1_yr_ago_GtC - C_in_ocean_1_yr_ago_GtC_LV2) / C_in_ocean_1_yr_ago_GtC_DL + (C_in_ocean_1_yr_ago_GtC_LV2 - C_in_ocean_1_yr_ago_GtC) / C_in_ocean_1_yr_ago_GtC_DL + (C_in_ocean_1_yr_ago_GtC - C_in_ocean_1_yr_ago_GtC_LV2) / C_in_ocean_1_yr_ago_GtC_DL),
            D(C_in_ocean_1_yr_ago_GtC_LV2) ~ - ((C_in_ocean_1_yr_ago_GtC_LV2 - C_in_ocean_1_yr_ago_GtC_LV1) / C_in_ocean_1_yr_ago_GtC_DL + ((C_in_ocean_1_yr_ago_GtC_LV2 - C_in_ocean_1_yr_ago_GtC_LV1) / C_in_ocean_1_yr_ago_GtC_DL + (C_in_ocean_1_yr_ago_GtC_LV1 - C_in_ocean_1_yr_ago_GtC_LV2) / C_in_ocean_1_yr_ago_GtC_DL)),
            D(C_in_ocean_1_yr_ago_GtC_LV1) ~ -((C_in_ocean_1_yr_ago_GtC_LV1 - Total_carbon_in_ocean_GtC) / C_in_ocean_1_yr_ago_GtC_DL + ((C_in_ocean_1_yr_ago_GtC_LV1 - Total_carbon_in_ocean_GtC) / C_in_ocean_1_yr_ago_GtC_DL + (Total_carbon_in_ocean_GtC - C_in_ocean_1_yr_ago_GtC_LV1) / C_in_ocean_1_yr_ago_GtC_DL)),
            0 ~ C_in_ocean_1_yr_ago_GtC_DL - 0.3333333333333333,
            0 ~ Atmos_heat_used_for_melting_last_year_1_yr - Atmos_heat_used_for_melting_last_year_1_yr_LV,
            D(Atmos_heat_used_for_melting_last_year_1_yr_LV) ~ Atmos_heat_used_for_melting_1_yr - Atmos_heat_used_for_melting_last_year_1_yr,
            0 ~ Ocean_heat_used_for_melting_last_year_ZJ_yr - Ocean_heat_used_for_melting_last_year_ZJ_yr_LV,
            D(Ocean_heat_used_for_melting_last_year_ZJ_yr_LV) ~ Ocean_heat_used_for_melting_ZJ_yr - Ocean_heat_used_for_melting_last_year_ZJ_yr,
            0 ~ C_in_atm_1_yr_ago_GtC + -C_in_atm_1_yr_ago_GtC_LV3 / C_in_atm_1_yr_ago_GtC_DL,
            D(C_in_atm_1_yr_ago_GtC_LV3) ~ C_in_atm_1_yr_ago_GtC_RT2 - C_in_atm_1_yr_ago_GtC,
            0 ~ C_in_atm_1_yr_ago_GtC_RT2 + -C_in_atm_1_yr_ago_GtC_LV2 / C_in_atm_1_yr_ago_GtC_DL,
            D(C_in_atm_1_yr_ago_GtC_LV2) ~ C_in_atm_1_yr_ago_GtC_RT1 - C_in_atm_1_yr_ago_GtC_RT2,
            0 ~ C_in_atm_1_yr_ago_GtC_RT1 + -C_in_atm_1_yr_ago_GtC_LV1 / C_in_atm_1_yr_ago_GtC_DL,
            D(C_in_atm_1_yr_ago_GtC_LV1) ~ C_in_atmosphere_GtC - C_in_atm_1_yr_ago_GtC_RT1,
            0 ~ C_in_atm_1_yr_ago_GtC_DL - 0.3333333333333333,
            0 ~
              All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC +
              -All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_LV3 / All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_DL,
            D(All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_LV3) ~
              All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_RT2 - All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC,
            0 ~
              All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_RT2 +
              -All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_LV2 / All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_DL,
            D(All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_LV2) ~
              All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_RT1 - All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_RT2,
            0 ~
              All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_RT1 +
              -All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_LV1 / All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_DL,
            D(All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_LV1) ~ All_C_taken_out_due_to_change_in_land_use_GtC - All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_RT1,
            0 ~ All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_DL - 0.3333333333333333,
            0 ~ (aux_1____Temp_gradient_minus_1___slope_ - 1.0) - Temp_gradient_minus_1___slope,
            0 ~
              Actual_time_to_degrade_all_GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_yr +
              -GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass / Sum_outflows_GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass_yr,
            0 ~
              Actual_time_to_degrade_all_NF_Dead_biomass__litter_and_soil_organic_matter_SOM_yr +
              -NF_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass / NF_Sum_outflows_GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass_yr,
            0 ~
              Actual_time_to_degrade_all_TROP_Dead_biomass__litter_and_soil_organic_matter_SOM_yr +
              -TROP_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass / TROP_Sum_outflows_GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass_yr,
            0 ~
              Actual_time_to_degrade_all_TUNDRA_Dead_biomass__litter_and_soil_organic_matter_SOM_yr - ESCIMO_ZIDZ(
                TUNDRA_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass,
                TUNDRA_Sum_outflows_GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass_yr,
              ),
            0 ~ Aerosol_anthropogenic_emissions - ESCIMO_IF_THEN_ELSE(Time >= 2015, 0.225Future_shape_of_anthropogenic_aerosol_emissions, Historical_aerosol_emissions_anthro),
            0 ~ Albedo_Antartic - ESCIMO_IF_THEN_ELSE(Time > 2020, 0.7, 0.7),
            0 ~ Albedo_glacier - ESCIMO_IF_THEN_ELSE(Time > 2020, 0.4, 0.4),
            0 ~
              (
                (
                  (
                    (
                      ((Albedo_land_biomes + (-0.24DESERT_Mkm2) / Area_of_land_Mkm2 + (-Albedo_URBAN * Urban_Mkm2) / Area_of_land_Mkm2) - Contrib_of_BARREN_land_to_albedo_land) -
                      Contrib_of_GRASS_to_albedo_land
                    ) - Contrib_of_ICE_ON_LAND_to_albedo_land
                  ) - Contrib_of_NF_to_albedo_land
                ) - Contrib_of_TROP_to_albedo_land
              ) - Contrib_of_TUNDRA_to_albedo_land,
            0 ~ (Albedo_ocean_with_arctic_ice_changes - 0.7Arctic_as_fraction_of_ocean) - 0.065Open_water_as_frac_of_ocean_area,
            0 ~ Albedo_URBAN - 0.15,
            0 ~
              All_C_taken_out_due_to_change_in_land_use_GtC -
              0.5 * (GRASS_land_taken_out_of_use_GtBiomass + NF_land_taken_out_of_use_GtBiomass + TROP_land_taken_out_of_use_GtBiomass + TUNDRA_land_taken_out_of_use_GtBiomass),
            0 ~
              (
                ((All_CH4_emissions_GtC_yr - CH4_release___capture_from_permafrost_area_loss___gain_GtC_yr) - Emissions_of_CH4_1850_to_2100_GtC_yr_with_IPCC_Fig_pg_1037_Exp) -
                Methanehydrate_experimental_release_GtC__yr
              ) - Natural_CH4_emissions,
          ]
        end
        function generateEquations5()
          println("#Equation generated:" * "50" * "in: " * "generateEquations5")
          [
            0 ~ (ALL_clouds_net_effect__pos_warming__neg_cooling__W_m2 - HI_clouds_net_effect__pos_warming__neg_cooling__W_m2) - LO_clouds_net_effect__pos_warming__neg_cooling__W_m2,
            0 ~ All_Human_activity_emissions_GtCO2e_yr_Base_for_tipping_point_search - var"combi_All_Human_activity_emissions_GtCO2e_yr_Base_for_tipping_point_search_y[1]",
            0 ~
              (
                (((All_Human_activity_emissions_GtCO2e_yr - Human_activity_CH4_emissions_GtCO2e_yr) - Kyoto_Flour_emissions_GtCO2e_yr) - Man_made_fossil_C_emissions_GtCO2e_yr) -
                Montreal_emissions_GtCO2e_yr
              ) - N2O_man_made_emissions_GtCO2e_yr,
            0 ~ (All_N2O_emissions_MtN2O_yr - 9.0) - N2O_man_made_emissions_exp_12a,
            0 ~ Annual_flux_of_C_to_biomass_GtC_pr_yr - Net_C_flow_from_atm_to_biomass_GtC_pr_yr,
            0 ~ Annual_glacial_ice_losing__pos__or_gaining__neg__GtIce_yr - 0.9167Glacial_ice_melting__pos__or_freezing__neg__km3_yr,
            0 ~ Annual_release_of_C_from_permafrost_GtC_y - CH4_release___capture_from_permafrost_area_loss___gain_GtC_yr,
            0 ~ Antarctic_ice_area_decrease_Mkm2_pr_yr + (-1.0e-6Antarctic_ice_melting_km3_yr) / Avg_thickness_Antarctic_km,
            0 ~ Antarctic_ice_area_increase_Mkm2_pr_yr + (-1.0e-6Antarctic_ice_freezing_km3_yr) / Avg_thickness_Antarctic_km,
            0 ~ Antarctic_ice_area_km2 + -Antarctic_ice_volume_km3 / Avg_thickness_Antarctic_km,
            0 ~
              Antarctic_ice_freezing_km3_yr -
              ESCIMO_IF_THEN_ELSE(Antarctic_ice_melting__pos__or_freezing__neg__km3_yr < 0.0, -Antarctic_ice_melting__pos__or_freezing__neg__km3_yr, 0.0),
            0 ~ Antarctic_ice_losing__pos__or_gaining__neg__GtIce_yr - 0.9167Antarctic_ice_melting__pos__or_freezing__neg__km3_yr,
            0 ~
              Antarctic_ice_melting__pos__or_freezing__neg__km3_yr +
              (
                -Antarctic_ice_volume_km3 *
                Effect_of_heat_in_atm_on_melting_ice__cut_off_ *
                Effect_of_temp_on_melting_antarctic_ice *
                Melting_constraint_from_the_heat_in__ocean__surface_reservoir *
                Melting_constraint_from_the_heat_in_atmosphere_reservoir_fraction *
                Snowball_earth_cutoff
              ) / Effective_time_to_melt_or_freeze_antarctic_ice_at_the_reference_delta_temp,
            0 ~ Antarctic_ice_melting_as_water_km3_yr - 0.916Antarctic_ice_melting__pos__or_freezing__neg__km3_yr,
            0 ~
              Antarctic_ice_melting_km3_yr -
              ESCIMO_IF_THEN_ELSE(Antarctic_ice_melting__pos__or_freezing__neg__km3_yr > 0.0, Antarctic_ice_melting__pos__or_freezing__neg__km3_yr, 0.0),
            0 ~ Anthropogenic_aerosol_forcing + 1.325 * Emissions_of_aerosols_1850_to_2100_with_IPCC_Fig_pg_1037_Exp,
            0 ~ Arctic_as_fraction_of_ocean + -Arctic_ice__on_sea__area_km2 / Ocean_area_km2,
            0 ~ Arctic_freezing_cutoff - var"combi_Arctic_freezing_cutoff_y[1]",
            0 ~ (Arctic_ice_area_max_km2 + Land_area_km2) - 5.1e8,
            0 ~ Arctic_ice_area_Mkm2 - 1.0e-6Arctic_ice__on_sea__area_km2,
            0 ~
              Arctic_ice_melting__pos__or_freezing__neg__km2_yr +
              (
                -Arctic_freezing_cutoff *
                Arctic_ice__on_sea__area_km2 *
                Effect_of_temp_on_melting_or_freezing_of_Arctic_ice *
                Melting_constraint_from_the_heat_in__ocean__surface_reservoir *
                Melting_constraint_from_the_heat_in_atmosphere_reservoir_fraction
              ) / Effective_time_to_melt_Arctic_ice_at_the_reference_delta_temp,
            0 ~
              Area_covered_by_high_clouds -
              (0.2 * (1.0 + Sensitivity_of_high_cloud_coverage_to_temp * (Temp_surface_current_divided_by_value_in_1850_K_K - 1.0))) *
              ESCIMO_IF_THEN_ELSE(Time > 2015, 1.0, 1.0),
            0 ~ Area_covered_by_low_clouds - (0.4 * (1.0 + 58.0 * (Temp_surface_current_divided_by_value_in_1850_K_K - 1.0))) * ESCIMO_IF_THEN_ELSE(Time > 2015, 1.0, 1.0),
            0 ~ Area_equivalent_of_linear_retreat_km2_yr - 12425.0,
            0 ~ Area_of_earth_Mkm2 - 510.0,
            0 ~ Area_of_land_Mkm2 - 0.30000000000000004Area_of_earth_Mkm2,
            0 ~ Atmos_heat_used_for_melting_1_yr + -Heat_withdrawn_from_atm_by_melting__pos____added__neg__by_freezing_ice_ZJ_yr / Heat_in_atmosphere_ZJ,
            0 ~ Avg_C_concentration_in_top_layer + -Carbon_in_top_ocean_layer_GtC / (Cold_surface_water_volume_Gm3 + Warm_surface_water_volume_Gm3),
            0 ~ Avg_CC_in_ocean_top_layer_ymoles_per_litre - Avg_C_concentration_in_top_layer * UNIT_converter_GtC_Gm3_to_ymoles_litre,
            0 ~ Avg_CO2_conc_in_ocean_top_layer_in_ppm - 0.127044Avg_CC_in_ocean_top_layer_ymoles_per_litre,
            0 ~ (Avg_earths_surface_albedo - 0.30000000000000004Albedo_land_biomes) - 0.7Albedo_ocean_with_arctic_ice_changes,
            0 ~ Avg_thickness_Antarctic_km - ESCIMO_IF_THEN_ELSE(Time > 2020, 2.14, 2.14),
            0 ~ Avg_thickness_glacier_km - ESCIMO_IF_THEN_ELSE(Time > 2020, 0.23, 0.23),
            0 ~ Avg_volcanic_activity_GtC_yr - 2.8Volcanic_aerosols_emissions,
            0 ~ (Barren_land_Mkm2 + Sum_biomes_Mkm2 + Urban_Mkm2) - 0.30000000000000004Area_of_earth_Mkm2,
            0 ~ (BARREN_land_normal_albedo_Mkm2 + BARREN_land_white_Mkm2) - Barren_land_Mkm2,
            0 ~ BARREN_land_white_Mkm2 - ESCIMO_IF_THEN_ELSE(false, EXP_12e_white_surfaces_ease_in, 0.0),
            0 ~ BB_radiation_at_atm_temp_in_atm_W_m2 + -BB_radiation_at_Temp_in_atm_ZJ_yr / UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~ BB_radiation_at_surface_temp_ZJ_yr - (5.67037e-8UNIT_conversion_W_m2_earth_to_ZJ_yr) * Temp_surface_average_K^4.0,
            0 ~ BB_radiation_at_Temp_in_atm_ZJ_yr - (5.67037e-8UNIT_conversion_W_m2_earth_to_ZJ_yr) * Temp_atm_average_K^4.0,
            0 ~ Blocked_by_CH4 - Fraction_blocked_by_CH4_spectrum,
            0 ~ Blocked_by_CO2 - Fraction_blocked_by_CO2_spectrum,
            0 ~
              Blocked_by_H20 - ESCIMO_IF_THEN_ELSE(
                false,
                Blocked_by_H20_Table_lookup,
                ESCIMO_IF_THEN_ELSE(false, Blocked_by_H2O_hist_and_fut, Blocked_by_h20_poly_used),
              ),
            0 ~ (Blocked_by_H20_future_linear_equ + Intercept_blocked_by_H20_future_equ) - Humidity_of_atmosphere_current_g_kg * Slope_blocked_by_H20_future_equ,
            0 ~
              ((Blocked_by_H20_future_poly_equ + Humidity_of_atmosphere_current_g_kg * exp1 + exp3 * Humidity_of_atmosphere_current_g_kg^3.0) - exp0) - exp2 * Humidity_of_atmosphere_current_g_kg^2.0,
            0 ~
              (((Blocked_by_H20_future_poly_equ_dyn - exp0_dyn) - Humidity_of_atmosphere_current_g_kg * exp1_dyn) - exp2_dyn * Humidity_of_atmosphere_current_g_kg^2.0) -
              exp3_dyn * Humidity_of_atmosphere_current_g_kg^3.0,
            0 ~ (((Blocked_by_H20_future_poly_equ_dyn_0 - exp0_dyn) - 2.0 * exp1_dyn) - 4.0 * exp2_dyn) - 8.0 * exp3_dyn,
            0 ~ Blocked_by_H20_hist_Table_lookup - var"combi_Blocked_by_H20_hist_Table_lookup_y[1]",
            0 ~ Blocked_by_h20_poly_used - ESCIMO_IF_THEN_ELSE(false, Blocked_by_H2O_poly_dyn, Blocked_by_H2O_poly_equ),
            0 ~ Blocked_by_H20_Table_lookup - var"combi_Blocked_by_H20_Table_lookup_y[1]",
          ]
        end
        function generateEquations6()
          println("#Equation generated:" * "50" * "in: " * "generateEquations6")
          [
            0 ~ Blocked_by_H2O_hist_and_fut - ESCIMO_IF_THEN_ELSE(Time > 2020, Blocked_by_H20_future_linear_equ, Blocked_by_H20_hist_Table_lookup),
            0 ~ Blocked_by_H2O_poly_dyn - ESCIMO_IF_THEN_ELSE(Time > 2020, Blocked_by_H20_future_poly_equ_dyn, Blocked_by_H20_hist_Table_lookup),
            0 ~ Blocked_by_H2O_poly_equ - ESCIMO_IF_THEN_ELSE(Time > 2020, Blocked_by_H20_future_poly_equ, Blocked_by_H20_hist_Table_lookup),
            0 ~ Blocked_by_otherGHG - Fraction_blocked_by_other_GHG,
            0 ~ Blocking_multiplier_from_Kyoto_Flour - ifEq_tmp304,
            0 ~ Blocking_multiplier_from_Montreal_gases - ifEq_tmp305,
            0 ~ (Blocking_multiplier_from_N2O - 1.0) - 0.1 * (N2O_concentration_ppb / Model_N2O_concentration_in_1850_ppb - 1.0),
            0 ~ (Blocking_of_LW_rad_by_clouds - LW_HI_cloud_radiation) - LW_LO_cloud_radiation,
            0 ~ C_absorption_by_ocean_from_atm_for_accumulation - C_diffusion_into_ocean_from_atm,
            0 ~
              C_diffusion_into_ocean_from_atm -
              Guldberg_Waage_air_sea_formulation * NatEvent_d__slowing_down_ocean_circulation_from_2015 * Ocean_circulation_slowdown_from_Greenland_ice_sliding_into_the_Atlantic,
            0 ~ C_diffusion_into_ocean_from_atm_MtC_yr - 1000.0C_diffusion_into_ocean_from_atm,
            0 ~ (((C_in_biomass - C_in_GRASS_GtC) - C_in_NF_GtC) - C_in_TROP_GtC) - C_in_TUNDRA_GtC,
            0 ~ C_in_GRASS_DeadB_and_soil_GtC - 0.5GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass,
            0 ~ ((C_in_GRASS_GtC - 0.5GRASS_Biomass_locked_in_construction_material_GtBiomass) - C_in_GRASS_DeadB_and_soil_GtC) - C_in_GRASS_LB_GtC,
            0 ~ C_in_GRASS_LB_GtC - 0.5GRASS_Living_biomass_GtBiomass,
            0 ~ C_in_NF_DeadB_and_soil_GtC - 0.5NF_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass,
            0 ~ ((C_in_NF_GtC - C_in_NF_DeadB_and_soil_GtC) - C_in_NF_LB_GtC) - 0.5NF_Biomass_locked_in_construction_material_GtBiomass,
            0 ~ C_in_NF_LB_GtC - 0.5NF_Living_biomass_GtBiomass,
            0 ~ C_in_TROP_DeadB_and_soil_GtC - 0.5TROP_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass,
            0 ~ ((C_in_TROP_GtC - 0.5TROP_Biomass_locked_in_construction_material_GtBiomass) - C_in_TROP_DeadB_and_soil_GtC) - C_in_TROP_LB_GtC,
            0 ~ C_in_TROP_LB_GtC - 0.5TROP_Living_biomass_GtBiomass,
            0 ~ C_in_TUNDRA_DeadB_and_soil_GtC - 0.5TUNDRA_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass,
            0 ~ ((C_in_TUNDRA_GtC - 0.5TUNDRA_Biomass_locked_in_construction_material_GtBiomass) - C_in_TUNDRA_DeadB_and_soil_GtC) - C_in_TUNDRA_LB_GtC,
            0 ~ C_in_TUNDRA_LB_GtC - 0.5TUNDRA_Living_biomass_GtBiomass,
            0 ~
              C_release_from_permafrost_melting_as_CO2_GtC_yr -
              CH4_in_permafrost_area_melted___frozen_before_heat_constraint_GtC_yr *
              Melting_restraint_for_permafrost_from_heat_in_atmophere *
              SHUT_OFF_permafrost *
              (1.0 - Fraction_of_C_released_from_permafrost_released_as_CH4_dmnl),
            0 ~ (C_released_from_permafrost_as_either_CH4_or_CO2_in_GtC_yr - CH4_release___capture_from_permafrost_area_loss___gain_GtC_yr) - C_release_from_permafrost_melting_as_CO2_GtC_yr,
            0 ~ C_removal_rate_from_atm_for_nature_May_2020_GtC_y - ESCIMO_IF_THEN_ELSE(Time > 2020, 0.0, 0.0),
            0 ~ C_runoff_from_biomass_soil - 0.5 * (GRASS_runoff + NF_runoff + TROP_runoff + TUNDRA_runoff),
            0 ~ Carbon_captured_and_stored_GtC___yr - ESCIMO_IF_THEN_ELSE(false, min(EXP_12b_CCS_from_2015, Man_made_fossil_C_emissions_GtC_yr), 0.0),
            0 ~
              Carbon_concentration_in_cold_surface_ocean +
              -C_in_cold_surface_water_GtC / (Cold_surface_water_volume_Gm3 + Cumulative_ocean_volume_increase_due_to_ice_melting_km3 * Frac_vol_cold_ocean_0_to_100m_of_total),
            0 ~
              Carbon_concentration_in_CWTtB +
              -C_in_cold_water_trunk_downwelling_GtC / (Cold_water_volume_downwelling_Gm3 + Cumulative_ocean_volume_increase_due_to_ice_melting_km3 * Frac_vol_cold_ocean_downwelling_of_total),
            0 ~
              Carbon_concentration_in_deep_box_GtC_per_G_cubicM +
              -C_in_deep_water_volume_1km_to_bottom_GtC / (Deep_water_volume_1km_to_4km_Gm3 + Cumulative_ocean_volume_increase_due_to_ice_melting_km3 * Frac_vol_deep_ocean_of_total),
            0 ~
              Carbon_concentration_in_intermdiate_box_GtC_per_G_cubicM +
              -C_in_intermediate_upwelling_water_100m_to_1km_GtC /
              (Intermediate_upwelling_water_volume_100m_to_1km_Gm3 + Cumulative_ocean_volume_increase_due_to_ice_melting_km3 * Frac_vol_ocean_upwelling_of_total),
            0 ~
              Carbon_concentration_in_warm_surface +
              -C_in_warm_surface_water_GtC / (Warm_surface_water_volume_Gm3 + Cumulative_ocean_volume_increase_due_to_ice_melting_km3 * Frac_vol_warm_ocean_0_to_100m_of_total),
            0 ~ Carbon_flow_from_cold_surface_downwelling_Gtc_per_yr - 0.1534278858251045C_in_cold_surface_water_GtC,
            0 ~ Carbon_flow_from_cold_to_deep_GtC_per_yr + -C_in_cold_water_trunk_downwelling_GtC / Time_in_trunk,
            0 ~ Carbon_flow_from_deep - 0.0013515522577680467C_in_deep_water_volume_1km_to_bottom_GtC,
            0 ~ Carbon_flow_from_intermediate_to_surface_box_GtC_per_yr - 0.004730436098903958C_in_intermediate_upwelling_water_100m_to_1km_GtC,
            0 ~ Carbon_flow_from_warm_to_cold_surface_GtC_per_yr - 0.0381286460517787C_in_warm_surface_water_GtC,
            0 ~ Carbon_in_cold_ocean_0_to_100m_1850_GtC + (-2240.0Volume_cold_ocean_0_to_100m) / UNIT_converter_GtC_Gm3_to_ymoles_litre,
            0 ~ Carbon_in_cold_ocean_trunk_100m_to_bottom_1850_GtC + (-2240.0Volume_cold_ocean_downwelling_100m_to_bottom) / UNIT_converter_GtC_Gm3_to_ymoles_litre,
            0 ~ Carbon_in_ocean_deep_1k_to_bottom_ocean_1850_GtC + (-2240.0Volume_ocean_deep_1km_to_bottom) / UNIT_converter_GtC_Gm3_to_ymoles_litre,
            0 ~ Carbon_in_ocean_upwelling_100m_to_1km_1850_GtC + (-2240.0Volume_ocean_upwelling_100m_to_1km) / UNIT_converter_GtC_Gm3_to_ymoles_litre,
            0 ~ (Carbon_in_top_ocean_layer_1850_GtC - Carbon_in_cold_ocean_0_to_100m_1850_GtC) - Carbon_in_warm_ocean_0_to_100m_1850_GtC,
            0 ~ (Carbon_in_top_ocean_layer_GtC - C_in_cold_surface_water_GtC) - C_in_warm_surface_water_GtC,
            0 ~ Carbon_in_warm_ocean_0_to_100m_1850_GtC + (-2240.0Volume_warm_ocean_0_to_100m) / UNIT_converter_GtC_Gm3_to_ymoles_litre,
            0 ~ CC_in_cold_downwelling_ymoles_per_litre - Carbon_concentration_in_CWTtB * UNIT_converter_GtC_Gm3_to_ymoles_litre,
            0 ~ CC_in_cold_downwelling_ymoles_per_litre__dimensionless_ - CC_in_cold_downwelling_ymoles_per_litre,
            0 ~ CC_in_cold_surface_ymoles_per_litre - Carbon_concentration_in_cold_surface_ocean * UNIT_converter_GtC_Gm3_to_ymoles_litre,
            0 ~ CC_in_cold_surface_ymoles_per_litre__dimensionless_ - CC_in_cold_surface_ymoles_per_litre,
          ]
        end
        function generateEquations7()
          println("#Equation generated:" * "50" * "in: " * "generateEquations7")
          [
            0 ~ CC_in_deep_box_ymoles_per_litre - Carbon_concentration_in_deep_box_GtC_per_G_cubicM * UNIT_converter_GtC_Gm3_to_ymoles_litre,
            0 ~ CC_in_deep_box_ymoles_per_litre__dimensionless_ - CC_in_deep_box_ymoles_per_litre,
            0 ~ CC_in_intermediate_box_ymoles_per_litre - Carbon_concentration_in_intermdiate_box_GtC_per_G_cubicM * UNIT_converter_GtC_Gm3_to_ymoles_litre,
            0 ~ CC_in_intermediate_box_ymoles_per_litre__dimensionless_ - CC_in_intermediate_box_ymoles_per_litre,
            0 ~ CC_in_warm_surface_ymoles_per_litre - Carbon_concentration_in_warm_surface * UNIT_converter_GtC_Gm3_to_ymoles_litre,
            0 ~ CC_in_warm_surface_ymoles_per_litre__dimensionless_ - CC_in_warm_surface_ymoles_per_litre,
            0 ~
              (((CH4_all_emissions_GtC_yr - CH4_release___capture_from_permafrost_area_loss___gain_GtC_yr) - Human_activity_CH4_emissions) - Methanehydrate_experimental_release_GtC__yr) -
              Natural_CH4_emissions,
            0 ~ CH4_concentration_ppb - MODEL_CH4_in_atm_in_ppb,
            0 ~ CH4_conversion_to_CO2_and_H2O - 0.136986301369863C_in_atmosphere_in_form_of_CH4,
            0 ~ CH4_emissions_before_co2e_exp - ESCIMO_IF_THEN_ELSE(false, 0.2, Emissions_of_CH4_1850_to_2100_GtC_yr_with_IPCC_Fig_pg_1037_Exp),
            0 ~
              CH4_emissions_CO2e_after_exp - ESCIMO_IF_THEN_ELSE(
                true,
                CH4_emissions_before_co2e_exp,
                ESCIMO_IF_THEN_ELSE(
                  false,
                  CH4_emissions_from_CO2e_C_Roads,
                  ESCIMO_IF_THEN_ELSE(
                    false,
                    CH4_emissions_from_CO2e_CAT,
                    ESCIMO_IF_THEN_ELSE(
                      false,
                      (0.01CH4_emissions_pct_contribution_to_Total_CO2e) * Total_CO2e_emissions_as_f_peak__GtCO2e_yr * UNIT_conversion_for_CH4_from_CO2e_to_C,
                      CH4_emissions_before_co2e_exp,
                    ),
                  ),
                ),
              ),
            0 ~ CH4_emissions_CO2e_after_exp_12a - ESCIMO_IF_THEN_ELSE(true, CH4_emissions_CO2e_after_exp, CH4_emissions_CO2e_after_exp * Exp_12a_reduction_in_emissions),
            0 ~ CH4_emissions_from_wetlands_destruction - CH4_per_sqkm_of_wetlands * Rate_of_destruction_of_wetlands,
            0 ~
              CH4_in_permafrost_area_melted___frozen_before_heat_constraint_GtC_yr -
              (4.8e-5Area_equivalent_of_linear_retreat_km2_yr) *
              Effect_of_temp_on_permafrost_melting_dmnl *
              Permafrost_melting_cutoff *
              Slowing_of_recapture_of_CH4_dmnl *
              ESCIMO_IF_THEN_ELSE(Time > 2015, 1.0, 1.0),
            0 ~ CH4_in_the_atmosphere_converted_to_CO2 - CH4_conversion_to_CO2_and_H2O,
            0 ~ CH4_per_sqkm_of_wetlands - 1.2e-5,
            0 ~
              CH4_release___capture_from_permafrost_area_loss___gain_GtC_yr -
              CH4_in_permafrost_area_melted___frozen_before_heat_constraint_GtC_yr *
              Fraction_of_C_released_from_permafrost_released_as_CH4_dmnl *
              Melting_restraint_for_permafrost_from_heat_in_atmophere *
              SHUT_OFF_permafrost,
            0 ~ (CO2_conc_atm_less_CO2_conc_sea + CO2_conc_in_cold_surface_water_in_ppm) - CO2_concentration_used__after_any_experiments__ppm,
            0 ~ CO2_conc_in_cold_surface_water_in_ppm - 0.127044CC_in_cold_surface_ymoles_per_litre,
            0 ~ CO2_conc_in_warm_surface_water_in_ppm - 0.127044CC_in_warm_surface_ymoles_per_litre,
            0 ~
              CO2_concentration_calculated_as_a_1pct_pa_exponential_increase_ppm -
              ESCIMO_IF_THEN_ELSE(Time > 20000000, CO2_ppm_value_at_When_to_sample * 1.01^Years_of_exponential_rise_dless, MODEL_CO2_concentration_in_atmosphere2_ppm),
            0 ~ CO2_concentration_ppm - CO2_concentration_used__after_any_experiments__ppm,
            0 ~
              CO2_concentration_used__after_any_experiments__ppm - ESCIMO_IF_THEN_ELSE(
                true,
                MODEL_CO2_concentration_in_atmosphere2_ppm,
                ESCIMO_IF_THEN_ELSE(false, dbl_CO2_exp, CO2_concentration_calculated_as_a_1pct_pa_exponential_increase_ppm),
              ),
            0 ~
              CO2_emissions_before_co2e_exp -
              RCPFossil_fuel_usage_cutoff *
              ESCIMO_IF_THEN_ELSE(false, Experimental_release_of_constant_fossil_C_emissions_GtC_yr, Emissions_of_CO2_1850_to_2100_GtC_yr_with_EXP_12a),
            0 ~
              CO2_emissions_CO2e_after_exp - ESCIMO_IF_THEN_ELSE(
                true,
                CO2_emissions_before_co2e_exp,
                ESCIMO_IF_THEN_ELSE(
                  false,
                  CO2_emissions_from_CO2e_C_Roads,
                  ESCIMO_IF_THEN_ELSE(
                    false,
                    CO2_emissions_from_CO2e_CAT,
                    ESCIMO_IF_THEN_ELSE(
                      false,
                      (0.01CO2_emissions_pct_contribution_to_Total_CO2e) * Total_CO2e_emissions_as_f_peak__GtCO2e_yr * UNIT_conversion_for_CO2_from_CO2e_to_C,
                      CO2_emissions_before_co2e_exp,
                    ),
                  ),
                ),
              ),
            0 ~
              CO2_flow_from_GRASS_to_atmosphere_GtC_yr -
              0.5 * (
                GRASS_Biomass_in_construction_material_being_burnt_GtBiomass_yr +
                GRASS_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr +
                GRASS_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr +
                GRASS_Dead_biomass_decomposing_GtBiomass_yr +
                GRASS_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr +
                GRASS_runoff +
                GRASS_soil_degradation_from_forest_fires_GtBiomass_yr
              ),
            0 ~
              CO2_flow_from_NF_to_atmosphere_GtC_yr -
              0.5 * (
                NF_Biomass_in_construction_material_being_burnt_GtBiomass_yr +
                NF_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr +
                NF_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr +
                NF_Dead_biomass_decomposing_GtBiomass_yr +
                NF_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr +
                NF_runoff +
                NF_soil_degradation_from_clear_cutting_GtBiomass_yr +
                NF_soil_degradation_from_forest_fires_GtBiomass_yr
              ),
            0 ~
              CO2_flow_from_TROP_to_atmosphere_GtC_yr -
              0.5 * (
                TROP_Biomass_in_construction_material_being_burnt_GtBiomass_yr +
                TROP_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr +
                TROP_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr +
                TROP_Dead_biomass_decomposing_GtBiomass_yr +
                TROP_NF_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr +
                TROP_runoff +
                TROP_soil_degradation_from_clear_cutting_GtBiomass_yr +
                TROP_soil_degradation_from_forest_fires_GtBiomass_yr
              ),
            0 ~
              CO2_flow_from_TUNDRA_to_atmosphere_GtC_yr -
              0.5 * (
                TUNDRA_Biomass_in_construction_material_being_burnt_GtBiomass_yr +
                TUNDRA_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr +
                TUNDRA_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr +
                TUNDRA_Dead_biomass_decomposing_GtBiomass_yr +
                TUNDRA_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr +
                TUNDRA_runoff +
                TUNDRA_soil_degradation_from_forest_fires_GtBiomass_yr
              ),
            0 ~ CO2_flux_from_atm_to_GRASS_for_new_growth_GtC_yr - 0.5GRASS_biomass_new_growing_GtBiomass___yr,
            0 ~ CO2_flux_from_atm_to_NF_for_new_growth_GtC_yr - 0.5NF_biomass_new_growing_GtBiomass___yr,
            0 ~ CO2_flux_from_atm_to_TROP_for_new_growth_GtC_yr - 0.5TROP_biomass_new_growing_GtBiomass___yr,
            0 ~ CO2_flux_from_atm_to_TUNDRA_for_new_growth_GtC_yr - 0.5TUNDRA_biomass_new_growing_GtBiomass___yr,
            0 ~ CO2_flux_GRASS_to_atm_Gtc_yr - CO2_flow_from_GRASS_to_atmosphere_GtC_yr,
            0 ~ CO2_flux_NF_to_atm_Gtc_yr - CO2_flow_from_NF_to_atmosphere_GtC_yr,
            0 ~ CO2_flux_TROP_to_atm_GtC_yr - CO2_flow_from_TROP_to_atmosphere_GtC_yr,
            0 ~ CO2_flux_TUNDRA_to_atm_Gtc_yr - CO2_flow_from_TUNDRA_to_atmosphere_GtC_yr,
            0 ~
              CO2_radiative_forcing_since_1850_using_Myhre_formula_W_pr_m2 -
              5.35 * ESCIMO_ln(CO2_concentration_used__after_any_experiments__ppm / CO2_concentration_in_1850_ppm),
            0 ~ Cold_dense_water_sinking_in_Sverdrup - (35.0NatEvent_d__slowing_down_ocean_circulation_from_2015) * Ocean_circulation_slowdown_from_Greenland_ice_sliding_into_the_Atlantic,
            0 ~ Concentration_of_C_in_ocean_top_layer_in_1850 + -Carbon_in_top_ocean_layer_1850_GtC / (Volume_cold_ocean_0_to_100m + Volume_warm_ocean_0_to_100m),
            0 ~ Contrib_of_BARREN_land_to_albedo_land + (-0.17BARREN_land_normal_albedo_Mkm2) / Area_of_land_Mkm2 + (-0.7BARREN_land_white_Mkm2) / Area_of_land_Mkm2,
            0 ~
              Contrib_of_GRASS_to_albedo_land +
              (-0.08GRASS_area_burnt_Mkm2) / Area_of_land_Mkm2 +
              (-0.3GRASS_deforested_Mkm2) / Area_of_land_Mkm2 +
              (-0.16 * ((GRASS_potential_area_Mkm2 - GRASS_area_burnt_Mkm2) - GRASS_deforested_Mkm2)) / Area_of_land_Mkm2,
            0 ~
              Contrib_of_ICE_ON_LAND_to_albedo_land +
              (-7.0e-7Greenland_ice_area_km2) / Area_of_land_Mkm2 +
              ((-1.0e-6Albedo_Antartic) * Antarctic_ice_area_km2) / Area_of_land_Mkm2 +
              ((-1.0e-6Albedo_glacier) * Glacial_ice_area_km2) / Area_of_land_Mkm2,
            0 ~
              Contrib_of_NF_to_albedo_land +
              (-0.08 * (((NF_potential_area_Mkm2 - NF_area_burnt_Mkm2) - NF_area_clear_cut_Mkm2) - NF_area_deforested_Mkm2)) / Area_of_land_Mkm2 +
              (-0.13NF_area_burnt_Mkm2) / Area_of_land_Mkm2 +
              (-0.18NF_area_clear_cut_Mkm2) / Area_of_land_Mkm2 +
              (-0.18NF_area_deforested_Mkm2) / Area_of_land_Mkm2,
            0 ~
              Contrib_of_TROP_to_albedo_land +
              (-0.1TROP_area_burnt_Mkm2) / Area_of_land_Mkm2 +
              (-0.168TROP_area_clear_cut_Mkm2) / Area_of_land_Mkm2 +
              (-0.168TROP_area_deforested_Mkm2) / Area_of_land_Mkm2 +
              (-0.14 * ((TROP_potential_area_Mkm2 - TROP_area_burnt_Mkm2) - TROP_area_deforested_Mkm2)) / Area_of_land_Mkm2,
            0 ~
              Contrib_of_TUNDRA_to_albedo_land +
              (-0.23TUNDRA_area_burnt_Mkm2) / Area_of_land_Mkm2 +
              (-0.23TUNDRA_deforested_Mkm2) / Area_of_land_Mkm2 +
              (-0.23 * ((Tundra_potential_area_Mkm2 - TUNDRA_area_burnt_Mkm2) - TUNDRA_deforested_Mkm2)) / Area_of_land_Mkm2,
            0 ~ Contribution_to_forcing_by_CH4 + -Blocked_by_CH4 / Frac_blocked_by_ALL_GHG,
            0 ~ Contribution_to_forcing_by_CO2 + -Blocked_by_CO2 / Frac_blocked_by_ALL_GHG,
            0 ~ Contribution_to_forcing_by_H2O + -Blocked_by_H20 / Frac_blocked_by_ALL_GHG,
            0 ~ Contribution_to_forcing_by_othGHG + -Blocked_by_otherGHG / Frac_blocked_by_ALL_GHG,
          ]
        end
        function generateEquations8()
          println("#Equation generated:" * "50" * "in: " * "generateEquations8")
          [
            0 ~ Convection_aka_sensible_heat_flow - Convection_as_f_of_temp_ZJ_yr,
            0 ~ Convection_aka_sensible_heat_flow_W_m2 + -Convection_aka_sensible_heat_flow / UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~ Convection_as_f_of_temp_ZJ_yr - (0.071Incoming_solar_in_1850_ZJ_yr) * (1.0 + 2.5 * (Temp_surface_current_divided_by_value_in_1850_K_K - 1.0)),
            0 ~ Conversion_constant_GtC_to_ppm + -600.0 / CO2_concentration_in_1850_ppm,
            0 ~ Conversion_constant_heat_ocean_deep_to_temp - 5.119803399549458e-7Temp__ocean__deep_in_1850_in_K,
            0 ~ Conversion_heat_atm_to_temp - 0.2674446946873751,
            0 ~ Conversion_heat_surface_to_temp - 0.0114726,
            0 ~ dbl_CO2_exp - ESCIMO_IF_THEN_ELSE(Time > 20000000, 2.0CO2_ppm_value_at_When_to_sample, MODEL_CO2_concentration_in_atmosphere2_ppm),
            0 ~ ((Deep_ocean__cold__volume - Cold_surface_water_volume_Gm3) - Cold_water_volume_downwelling_Gm3) - Deep_water_volume_1km_to_4km_Gm3,
            0 ~ (C_in_atmosphere_GtC_in_1850 + delta_C_in_atmosphere_GtC) - C_in_atmosphere_GtC,
            0 ~ (C_in_biomass_in_1850_GtC + delta_C_in_biomass_GtC) - C_in_biomass,
            0 ~ (Total_carbon_in_ocean_GtC_in_1850 + delta_C_in_ocean_GtC) - Total_carbon_in_ocean_GtC,
            0 ~ (CO2_concentration_in_1850_ppm + delta_CO2_concentration_since_1850_ppm) - CO2_concentration_used__after_any_experiments__ppm,
            0 ~ (Temp_ocean_deep_1850_degC + delta_Temp_deep_ocean_degC) - Temp__ocean__deep_in_C,
            0 ~ Depositing_of_C_to_sediment - 5.0e-5C_in_deep_water_volume_1km_to_bottom_GtC,
            0 ~ (Effect_of_acidification_on_NMPP - 1.0) - 5.0 * (ph_in_cold_downwelling_water / init_ph_in_cold_water - 1.0),
            0 ~ (Effect_of_C_concentration_on_NMPP - 1.0) - 1.4426950408889634 * ESCIMO_ln(Avg_C_concentration_in_top_layer / Concentration_of_C_in_ocean_top_layer_in_1850),
            0 ~ (Effect_of_CO2_on_new_biomass_growth - 1.0) - ESCIMO_ln(CO2_concentration_used__after_any_experiments__ppm / CO2_concentration_in_1850_ppm),
            0 ~ Effect_of_heat_in_atm_on_melting_ice__cut_off_ - var"combi_Effect_of_heat_in_atm_on_melting_ice__cut_off__y[1]",
            0 ~ (Effect_of_humidity_on_shifting_biomes - 1.0) - 5.0 * (Humidity_of_atmosphere_g_kg / Humidity_of_atmosphere_in_1850_g_kg - 1.0),
            0 ~
              Effect_of_population_and_urbanization_on_biomass_use -
              (0.1639344262295082Urbanzation_Effect_on_biomass_use) * ESCIMO_Population_Lookup_bn(var"combi_Population_Lookup_bn_y[1]"),
            0 ~ (Effect_of_temp_on_melting_antarctic_ice - 1.0) - 1.2 * (0.3333333333333333Temp_diff_relevant_for_melting_or_freezing_from_1850 - 1.0),
            0 ~ (Effect_of_temp_on_melting_greenland_ice - 1.0) - 0.1 * (Arctic_land_surface_temp_anomaly_compared_to_1850 - 1.0),
            0 ~ (Effect_of_temp_on_melting_greenland_ice_that_slid_into_the_ocean - 1.0) - 0.71 * (Temp_surface_C - 14.66500000000002),
            0 ~ (Effect_of_temp_on_melting_or_freezing_glacial_ice - 1.0) - Slope_temp_vs_glacial_ice_melting * (0.3333333333333333Temp_diff_relevant_for_melting_or_freezing_from_1850 - 1.0),
            0 ~ (Effect_of_temp_on_melting_or_freezing_of_Arctic_ice - 1.0) - 0.65 * (2.5Temp_diff_relevant_for_melting_or_freezing_arctic_ice_from_1850 - 1.0),
            0 ~ (Effect_of_temperature_on_fire_incidence_dimensionless - 1.0) - (0.1 * (0.07317965605561642Temp_surface_C - 1.0)) * ESCIMO_IF_THEN_ELSE(Time > 2015, 1.0, 1.0),
            0 ~ (Effect_of_temperature_on_new_biomass_growth_dimensionless + 0.5 * (0.07317965605561642Temp_surface_C - 1.0)) - 1.0,
            0 ~ (Effect_of_temperature_on_NMPP - 1.0) - 2.0 * (0.07317965605561642Temp_surface_C - 1.0),
            0 ~ Effective_time_to_melt_Arctic_ice_at_the_reference_delta_temp - 500.0 * ESCIMO_IF_THEN_ELSE(Time > 2015, 1.0, 1.0),
            0 ~ Effective_time_to_melt_glacial_ice_at_the_reference_delta_temp - 500.0 * ESCIMO_IF_THEN_ELSE(Time > 2015, 1.0, 1.0),
            0 ~ Effective_time_to_melt_greenland_ice_at_the_reference_delta_temp - 4000.0 * ESCIMO_IF_THEN_ELSE(Time > 2015, 1.0, 1.0),
            0 ~ Effective_time_to_melt_or_freeze_antarctic_ice_at_the_reference_delta_temp - 18000.0 * ESCIMO_IF_THEN_ELSE(Time > 2015, 1.0, 1.0),
            0 ~ Effective_Time_to_regrow_TROP_after_deforesting_yr + -10000.0 / TROP_deforestation_cutoff_effect,
            0 ~
              Emissions_of_aerosols_1850_to_2100_with_IPCC_Fig_pg_1037_Exp - ESCIMO_IF_THEN_ELSE(
                Time <= 2010,
                Aerosol_anthropogenic_emissions,
                ESCIMO_IF_THEN_ELSE(
                  true,
                  Aerosol_anthropogenic_emissions,
                  ESCIMO_IF_THEN_ELSE(false, 0.0, ESCIMO_IF_THEN_ELSE(false, Aerosol_anthropogenic_emissions_in_2010, 0.0)),
                ),
              ),
            0 ~
              Emissions_of_anthro_CH4_1850_to_2100_GtC_yr - ESCIMO_IF_THEN_ELSE(
                true,
                Emissions_of_anthro_CH4_1850_to_2100_RCP_with_JR_shape_forecast_GtC_yr,
                ESCIMO_IF_THEN_ELSE(
                  false,
                  Emissions_of_anthro_CH4_1850_to_2100_RCP_with_Xpct_annual_incr_forecast_GtC_yr,
                  ESCIMO_IF_THEN_ELSE(
                    false,
                    Emissions_of_anthro_CH4_1850_to_2100_linearly_reduced_from_2015_GtC_yr,
                    ESCIMO_IF_THEN_ELSE(
                      false,
                      Emissions_of_anthro_CH4_1850_to_2100_RCP3_GtC_yr,
                      ESCIMO_IF_THEN_ELSE(
                        false,
                        Emissions_of_anthro_CH4_1850_to_2100_RCP45_GtC_yr,
                        ESCIMO_IF_THEN_ELSE(
                          false,
                          Emissions_of_anthro_CH4_1850_to_2100_RCP6_GtC_yr,
                          ESCIMO_IF_THEN_ELSE(false, Emissions_of_anthro_CH4_1850_to_2100_RCP85_GtC_yr, 0.0),
                        ),
                      ),
                    ),
                  ),
                ),
              ),
            0 ~
              Emissions_of_anthro_CH4_1850_to_2100_linearly_reduced_from_2015_GtC_yr - ESCIMO_IF_THEN_ELSE(
                Time >= 2015,
                0.303 * min(1.0, max(0.0, (2050.0 - Time) / Years_still_needed_to_reach_zero_emission_goal_yr)),
                Emissions_of_anthro_CH4_1850_to_2100_RCP_with_JR_shape_forecast_GtC_yr,
              ),
            0 ~
              Emissions_of_anthro_CH4_1850_to_2100_RCP_with_JR_shape_forecast_GtC_yr -
              Emissions_of_anthro_CH4_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_MtCH4_yr * UNIT_conversion_from_MtCH4_to_GtC,
            0 ~
              Emissions_of_anthro_CH4_1850_to_2100_RCP_with_Xpct_annual_incr_forecast_GtC_yr -
              ESCIMO_IF_THEN_ELSE(Time >= 2015, 0.303, Emissions_of_anthro_CH4_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_MtCH4_yr * UNIT_conversion_from_MtCH4_to_GtC),
            0 ~ Emissions_of_anthro_CH4_1850_to_2100_RCP3_GtC_yr - Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP3_MtCH4_yr * UNIT_conversion_from_MtCH4_to_GtC,
            0 ~ Emissions_of_anthro_CH4_1850_to_2100_RCP45_GtC_yr - Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP45_MtCH4_yr * UNIT_conversion_from_MtCH4_to_GtC,
            0 ~ Emissions_of_anthro_CH4_1850_to_2100_RCP6_GtC_yr - Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP6_MtCH4_yr * UNIT_conversion_from_MtCH4_to_GtC,
            0 ~ Emissions_of_anthro_CH4_1850_to_2100_RCP85_GtC_yr - Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP85_MtCH4_yr * UNIT_conversion_from_MtCH4_to_GtC,
            0 ~
              Emissions_of_anthro_CO2_1850_to_2100_linearly_reduced_from_2015_GtC_yr - ESCIMO_IF_THEN_ELSE(
                Time >= 2015,
                10.0 * min(1.0, max(0.0, (2050.0 - Time) / Years_still_needed_to_reach_zero_emission_goal_yr)),
                Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_GtC_yr,
              ),
            0 ~
              Emissions_of_anthro_CO2_1850_to_2100_RCP_with_Xpct_annual_incr_forecast_GtC_yr -
              ESCIMO_IF_THEN_ELSE(Time >= 2015, 10.0, Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_GtC_yr),
            0 ~
              Emissions_of_CH4_1850_to_2100_GtC_yr_with_IPCC_Fig_pg_1037_Exp - ESCIMO_IF_THEN_ELSE(
                Time <= 2010,
                Emissions_of_anthro_CH4_1850_to_2100_GtC_yr,
                ESCIMO_IF_THEN_ELSE(
                  true,
                  Emissions_of_anthro_CH4_1850_to_2100_GtC_yr,
                  ESCIMO_IF_THEN_ELSE(false, 0.0, ESCIMO_IF_THEN_ELSE(false, CO4_emissions_in_2010, 0.0)),
                ),
              ),
            0 ~
              Emissions_of_CO2_1850_to_2100_GtC_yr_with_EXP_12a - ESCIMO_IF_THEN_ELSE(
                true,
                Emissions_of_CO2_1850_to_2100_GtC_yr_with_IPCC_Fig_pg_1037_Exp,
                Emissions_of_CO2_1850_to_2100_GtC_yr_with_IPCC_Fig_pg_1037_Exp * Exp_12a_reduction_in_emissions,
              ),
            0 ~
              Emissions_of_CO2_1850_to_2100_GtC_yr_with_IPCC_Fig_pg_1037_Exp - ESCIMO_IF_THEN_ELSE(
                Time <= 2010,
                Emissions_of_CO2_1850_to_2100_GtC_yr,
                ESCIMO_IF_THEN_ELSE(
                  true,
                  Emissions_of_CO2_1850_to_2100_GtC_yr,
                  ESCIMO_IF_THEN_ELSE(false, 0.0, ESCIMO_IF_THEN_ELSE(false, CO2_emissions_in_2010, 0.0)),
                ),
              ),
            0 ~
              Emissions_of_CO2_1850_to_2100_GtC_yr - ESCIMO_IF_THEN_ELSE(
                false,
                Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_GtC_yr,
                ESCIMO_IF_THEN_ELSE(
                  false,
                  Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_hist_with_JR__large_scale_CCS__exp_GtC_yr,
                  ESCIMO_IF_THEN_ELSE(
                    true,
                    Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_GtC_yr,
                    ESCIMO_IF_THEN_ELSE(
                      false,
                      Emissions_of_anthro_CO2_1850_to_2100_linearly_reduced_from_2015_GtC_yr,
                      ESCIMO_IF_THEN_ELSE(
                        false,
                        Emissions_of_anthro_CO2_1850_to_2100_RCP_with_Xpct_annual_incr_forecast_GtC_yr,
                        ESCIMO_IF_THEN_ELSE(
                          false,
                          Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP3_GtC_yr,
                          ESCIMO_IF_THEN_ELSE(
                            false,
                            Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP45_GtC_yr,
                            ESCIMO_IF_THEN_ELSE(
                              false,
                              Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP6_GtC_yr,
                              ESCIMO_IF_THEN_ELSE(false, Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP85_GtC_yr, 0.0),
                            ),
                          ),
                        ),
                      ),
                    ),
                  ),
                ),
              ),
            0 ~ Evaporation_aka_latent_heat_flow - Evaporation_as_f_of_temp_ZJ_yr,
          ]
        end
        function generateEquations9()
          println("#Equation generated:" * "50" * "in: " * "generateEquations9")
          [
            0 ~ Evaporation_aka_latent_heat_flow_W_m2 + -Evaporation_aka_latent_heat_flow / UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~ Evaporation_as_f_of_temp_ZJ_yr - (0.289Incoming_solar_in_1850_ZJ_yr) * (1.0 + 0.057999999999999996Temp_surface_anomaly_compared_to_1850_degC),
            0 ~
              Exogenous_sliding_of_Greenland_ice_into_the_ocean -
              ESCIMO_IF_THEN_ELSE(Greenland_ice_volume_on_Greenland_km3 < Greenland_slide_experiment_end_condition, 0.0, 1.0) *
              ESCIMO_IF_THEN_ELSE(Time > 3000000, 0.00510204081632653, 0.0),
            0 ~
              Exp_12a_reduction_in_emissions -
              ESCIMO_IF_THEN_ELSE(Time < 2015, 1.0, ESCIMO_IF_THEN_ELSE(Time > 2087, 0.0, Exp_12a_reduction_in_emissions_LOOKUP)),
            0 ~ Exp_12a_reduction_in_emissions_LOOKUP - var"combi_Exp_12a_reduction_in_emissions_LOOKUP_y[1]",
            0 ~ EXP_12b_CCS_from_2015 - var"combi_EXP_12b_CCS_from_2015_y[1]",
            0 ~ EXP_12c_stopping_TROP_deforestation_from_2015 - ESCIMO_IF_THEN_ELSE(Time > 2015, 0.0, 1.0),
            0 ~ EXP_12e_white_surfaces_ease_in - var"combi_EXP_12e_white_surfaces_ease_in_y[1]",
            0 ~ exp0 - 2.4523,
            0 ~ (exp0_dyn - 2.4523) - 0.343Slider_for_H2O_slope,
            0 ~ exp1 - 3.7148,
            0 ~ 3.7148 + exp1_dyn + 0.44Slider_for_H2O_slope,
            0 ~ exp2 - 1.8344,
            0 ~ (exp2_dyn - 1.8344) - 0.1794Slider_for_H2O_slope,
            0 ~ exp3 - 0.2842,
            0 ~ 0.2842 + exp3_dyn + 0.0227Slider_for_H2O_slope,
            0 ~ ((Experimental_doubling_of_constant_C_emissions + ESCIMO_STEP(t, 0.0, 30005.0)) - 1.0) - ESCIMO_STEP(t, 0.0, 30000.0),
            0 ~ Experimental_release_of_constant_fossil_C_emissions_GtC_yr - 4.0 * Experimental_doubling_of_constant_C_emissions,
            0 ~ (Experimental_release_of_methane + ESCIMO_STEP(t, 0.0, 2025.0)) - ESCIMO_STEP(t, 0.0, 2020.0),
            0 ~ f_M_1750_N_2010__for_ch4_forcing - 0.47 * ESCIMO_ln(1.0 + 3.045720946785617e-14 * N_2010^1.52 + 3.015e-5N_2010),
            0 ~ f_M_2010_N_cur_ - 0.47 * ESCIMO_ln(1.0 + (5.31e-15 * M_2010^2.52) * N_cur^1.52 + (1.5075e-5M_2010) * N_cur),
            0 ~ f_M_cur_N_2010_ - 0.47 * ESCIMO_ln(1.0 + (5.31e-15 * M_cur^2.52) * N_2010^1.52 + (1.5075e-5M_cur) * N_2010),
            0 ~ f_M2010_N_1750__for_n20_forcing - 0.47 * ESCIMO_ln(1.0 + 3.015e-5M_2010 + 1.5228604733928084e-14 * M_2010^2.52),
            0 ~
              (
                ((Flow_from_atm_to_biomass_GtC_pr_yr - CO2_flux_from_atm_to_GRASS_for_new_growth_GtC_yr) - CO2_flux_from_atm_to_NF_for_new_growth_GtC_yr) -
                CO2_flux_from_atm_to_TROP_for_new_growth_GtC_yr
              ) - CO2_flux_from_atm_to_TUNDRA_for_new_growth_GtC_yr,
            0 ~ (((Flow_from_biomass_to_atm_Gtc_pr_yr - CO2_flux_GRASS_to_atm_Gtc_yr) - CO2_flux_NF_to_atm_Gtc_yr) - CO2_flux_TROP_to_atm_GtC_yr) - CO2_flux_TUNDRA_to_atm_Gtc_yr,
            0 ~ Flow_of_cold_surface_water_welling_down_GcubicM_per_yr - Flow_of_cold_water_sinking_to_very_bottom_GcubicM_per_yr,
            0 ~ Flow_of_cold_water_sinking_to_very_bottom_GcubicM_per_yr - 31536.0Cold_dense_water_sinking_in_Sverdrup,
            0 ~ Flow_of_heat_to_atm_ZJ_yr - Net_heat_flow_to_atm_ZJ_yr__needed_for_comparisons_with_history_,
            0 ~ Flow_of_heat_to_deep_ocean - Net_heat_flow_ocean_from_surface_to_deep__ZJ_yr_,
            0 ~
              Flow_of_heat_to_deep_ocean_btw_72_and_08 -
              ESCIMO_IF_THEN_ELSE(Time >= 2008, 0.0, ESCIMO_IF_THEN_ELSE(Time >= 1972, Net_heat_flow_ocean_from_surface_to_deep__ZJ_yr_, 0.0)),
            0 ~ Flow_of_heat_to_surface_ocean - Net_flow_of_heat_into_surface,
            0 ~
              Flow_of_heat_to_surface_ocean_btw_1972_and_2008 -
              ESCIMO_IF_THEN_ELSE(Time >= 2008, 0.0, ESCIMO_IF_THEN_ELSE(Time >= 1972, Net_flow_of_heat_into_surface, 0.0)),
            0 ~ Flow_of_water_from_warm_to_cold_surface_Gcubicm_per_yr - Flow_of_cold_water_sinking_to_very_bottom_GcubicM_per_yr,
            0 ~ for_display_yr_on_yr_change_in_C_in_ocean_GtC_yr + yr_on_yr_change_in_C_in_ocean_GtC_yr,
            0 ~ Frac_atm_absorption - ESCIMO_IF_THEN_ELSE(Time > 2020, 0.220588, Hist_Frac_atm_absorption),
            0 ~ (((Frac_blocked_by_ALL_GHG - Blocked_by_CH4) - Blocked_by_CO2) - Blocked_by_H20) - Blocked_by_otherGHG,
            0 ~ ((Frac_blocked_by_ALL_GHG_LESS_watervapor - Blocked_by_CH4) - Blocked_by_CO2) - Blocked_by_otherGHG,
            0 ~ Frac_vol_cold_ocean_0_to_100m_of_total + -Volume_cold_ocean_0_to_100m / Volume_of_total_ocean_Gm3,
            0 ~ Frac_vol_cold_ocean_downwelling_of_total + -Volume_cold_ocean_downwelling_100m_to_bottom / Volume_of_total_ocean_Gm3,
            0 ~ Frac_vol_deep_ocean_of_total + -Volume_ocean_deep_1km_to_bottom / Volume_of_total_ocean_Gm3,
            0 ~ Frac_vol_ocean_upwelling_of_total + -Volume_ocean_upwelling_100m_to_1km / Volume_of_total_ocean_Gm3,
            0 ~ Frac_vol_warm_ocean_0_to_100m_of_total + -Volume_warm_ocean_0_to_100m / Volume_of_total_ocean_Gm3,
            0 ~ Fraction_blocked_by_CH4_spectrum - var"combi_Fraction_blocked_by_CH4_spectrum_y[1]",
            0 ~ Fraction_blocked_by_CO2_spectrum - var"combi_Fraction_blocked_by_CO2_spectrum_y[1]",
            0 ~ Fraction_blocked_by_other_GHG - 0.0398LW_Blocking_multiplier_from_other_GHG,
            0 ~ Fraction_GRASS_being_deforested_1_yr - GRASS_historical_deforestation_pct_yr,
            0 ~ Fraction_of_C_released_from_permafrost_released_as_CH4_dmnl - ESCIMO_IF_THEN_ELSE(Time < 2020, 1.0, 1.0),
            0 ~ Fraction_of_ocean_classified_as_cold_surface - 0.19999999999999996,
            0 ~ Fraction_TUNDRA_being_deforested_1_yr - TUNDRA_historical_deforestation_pct_yr,
            0 ~ Future_shape_of_anthropogenic_aerosol_emissions - var"combi_Future_shape_of_anthropogenic_aerosol_emissions_y[1]",
          ]
        end
        function generateEquations10()
          println("#Equation generated:" * "50" * "in: " * "generateEquations10")
          [
            0 ~ (Ga__BB_radiation_less_TOA_radiation_W_m2 + LW_TOA_radiation_from_atm_to_space_W_m2) - BB_radiation_at_atm_temp_in_atm_W_m2,
            0 ~ Glacial_ice_area_decrease_Mkm2_pr_yr + (-1.0e-6Glacial_ice_melting_km3_yr) / Avg_thickness_glacier_km,
            0 ~ Glacial_ice_area_increase_Mkm2_pr_yr + (-1.0e-6Glacial_ice_freezing_km3_yr) / Avg_thickness_glacier_km,
            0 ~ Glacial_ice_area_km2 + -Glacial_ice_volume_km3 / Avg_thickness_glacier_km,
            0 ~
              Glacial_ice_freezing_km3_yr -
              ESCIMO_IF_THEN_ELSE(Glacial_ice_melting__pos__or_freezing__neg__km3_yr < 0.0, -Glacial_ice_melting__pos__or_freezing__neg__km3_yr, 0.0),
            0 ~
              Glacial_ice_melting__pos__or_freezing__neg__km3_yr +
              (-Effect_of_heat_in_atm_on_melting_ice__cut_off_ * Effect_of_temp_on_melting_or_freezing_glacial_ice * Glacial_ice_volume_km3) /
              Effective_time_to_melt_glacial_ice_at_the_reference_delta_temp,
            0 ~ Glacial_ice_melting_as_water_km3_yr - 0.916Glacial_ice_melting__pos__or_freezing__neg__km3_yr,
            0 ~
              Glacial_ice_melting_km3_yr -
              ESCIMO_IF_THEN_ELSE(Glacial_ice_melting__pos__or_freezing__neg__km3_yr > 0.0, Glacial_ice_melting__pos__or_freezing__neg__km3_yr, 0.0),
            0 ~ GRASS_being_deforested_Mkm2_yr - Fraction_GRASS_being_deforested_1_yr * GRASS_with_normal_cover_Mkm2,
            0 ~ GRASS_being_harvested_Mkm2_yr + (-1000.0Use_of_GRASS_biomass_for_energy_GtBiomass_yr) / GRASS_living_biomass_densitiy_tBiomass_pr_km2,
            0 ~
              GRASS_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr -
              (0.001GRASS_living_biomass_densitiy_tBiomass_pr_km2) * (GRASS_being_deforested_Mkm2_yr + GRASS_being_harvested_Mkm2_yr + GRASS_burning_Mkm2_yr),
            0 ~ GRASS_Biomass_in_construction_material_being_burnt_GtBiomass_yr - 0.05GRASS_Biomass_locked_in_construction_material_GtBiomass,
            0 ~ GRASS_Biomass_in_construction_material_left_to_rot_GtBiomass_yr - 0.05GRASS_Biomass_locked_in_construction_material_GtBiomass,
            0 ~ GRASS_biomass_new_growing_GtBiomass___yr - 0.5GRASS_potential_less_actual_living_biomass_GtBiomass,
            0 ~ GRASS_burning_Mkm2_yr - (0.01 * Effect_of_temperature_on_fire_incidence_dimensionless) * GRASS_with_normal_cover_Mkm2,
            0 ~ GRASS_Dead_biomass_decomposing_GtBiomass_yr - 0.001GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass,
            0 ~ GRASS_DeadB_and_SOM_tB_per_km2 + (-1000.0GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass) / GRASS_with_normal_cover_Mkm2,
            0 ~ GRASS_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr - (0.001GRASS_DeadB_and_SOM_tB_per_km2) * GRASS_being_deforested_Mkm2_yr,
            0 ~ GRASS_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr - (0.0001GRASS_DeadB_and_SOM_tB_per_km2) * GRASS_being_harvested_Mkm2_yr,
            0 ~ GRASS_for_construction_use_GtBiomass_yr - Use_of_GRASS_biomass_for_construction_GtBiomass_yr,
            0 ~ GRASS_historical_deforestation_pct_yr - 0.001,
            0 ~ GRASS_land_taken_out_of_use_GtBiomass - (0.001GRASS_land_taken_out_of_use_Mkm2) * GRASS_living_biomass_densitiy_tBiomass_pr_km2,
            0 ~ ((GRASS_land_taken_out_of_use_Mkm2 - GRASS_area_burnt_Mkm2) - GRASS_area_harvested_Mkm2) - GRASS_deforested_Mkm2,
            0 ~ GRASS_living_biomass_densitiy_tBiomass_pr_km2 - (14500.0 * Effect_of_CO2_on_new_biomass_growth) * Effect_of_temperature_on_new_biomass_growth_dimensionless,
            0 ~ GRASS_Living_biomass_rotting_GtBiomass_yr - 0.01GRASS_Living_biomass_GtBiomass,
            0 ~ (GRASS_Living_biomass_GtBiomass + GRASS_potential_less_actual_living_biomass_GtBiomass) - GRASS_potential_living_biomass_GtBiomass,
            0 ~ GRASS_potential_living_biomass_GtBiomass - (0.001GRASS_living_biomass_densitiy_tBiomass_pr_km2) * (GRASS_potential_area_Mkm2 - GRASS_deforested_Mkm2),
            0 ~ GRASS_regrowing_after_being_burnt_Mkm2_yr - 0.1GRASS_area_burnt_Mkm2,
            0 ~ GRASS_regrowing_after_being_deforested_Mkm2_yr - 0.0125GRASS_deforested_Mkm2,
            0 ~ GRASS_regrowing_after_harvesting_Mkm2_yr - 0.1GRASS_area_harvested_Mkm2,
            0 ~ GRASS_runoff - 0.0005GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass,
            0 ~ GRASS_soil_degradation_from_forest_fires_GtBiomass_yr,
            0 ~ (GRASS_area_burnt_Mkm2 + GRASS_area_harvested_Mkm2 + GRASS_deforested_Mkm2 + GRASS_with_normal_cover_Mkm2) - GRASS_potential_area_Mkm2,
            0 ~ Greenland_ice_area_decrease_Mkm2_pr_yr - 7.407407407407406e-7Greenland_ice_melting_km3_yr,
            0 ~ Greenland_ice_area_increase_Mkm2_pr_yr - 7.407407407407406e-7Greenland_ice_freezing_km3_yr,
            0 ~ Greenland_ice_area_km2 - 0.7407407407407407Greenland_ice_volume_on_Greenland_km3,
            0 ~
              Greenland_ice_freezing_km3_yr -
              ESCIMO_IF_THEN_ELSE(Greenland_ice_melting__pos__or_freezing__neg__km3_yr < 0.0, -Greenland_ice_melting__pos__or_freezing__neg__km3_yr, 0.0),
            0 ~ Greenland_ice_losing__pos__or_gaining__neg__GtIce_yr - 0.9167Greenland_ice_melting__pos__or_freezing__neg__km3_yr,
            0 ~
              Greenland_ice_melting__pos__or_freezing__neg__km3_yr +
              (
                -Effect_of_heat_in_atm_on_melting_ice__cut_off_ *
                Effect_of_temp_on_melting_greenland_ice *
                Greenland_ice_volume_on_Greenland_km3 *
                Melting_constraint_from_the_heat_in_atmosphere_reservoir_fraction *
                Snowball_earth_cutoff
              ) / Effective_time_to_melt_greenland_ice_at_the_reference_delta_temp,
            0 ~ Greenland_ice_melting_as_water_km3_yr - 0.916Greenland_ice_melting__pos__or_freezing__neg__km3_yr,
            0 ~
              Greenland_ice_melting_km3_yr -
              ESCIMO_IF_THEN_ELSE(Greenland_ice_melting__pos__or_freezing__neg__km3_yr > 0.0, Greenland_ice_melting__pos__or_freezing__neg__km3_yr, 0.0),
            0 ~ Greenland_ice_melting_that_slid_into_the_ocean_km3_yr - Greenland_ice_sliding_into_the_ocean_km3_yr,
            0 ~ Greenland_ice_sliding_into_the_ocean_km3_yr - Exogenous_sliding_of_Greenland_ice_into_the_ocean * Greenland_ice_volume_on_Greenland_km3,
            0 ~
              Greenland_ice_that_slid_into_the_ocean_melting__pos__or_freezing__neg__km3_yr -
              (0.05 * Effect_of_temp_on_melting_greenland_ice_that_slid_into_the_ocean) *
              Greenland_ice_volume_that_slid_into_the_ocean_km3 *
              Melting_constraint_from_the_heat_in__ocean__surface_reservoir,
            0 ~ Guldberg_Waage_air_sea_formulation - (0.05555555555555555CO2_conc_atm_less_CO2_conc_sea) * Conversion_constant_GtC_to_ppm,
            0 ~ Heat_actually_gained___needed_for_freezing___unfreezing_of_permafrost_ZJ_yr - 0.35770833333333335CH4_release___capture_from_permafrost_area_loss___gain_GtC_yr,
            0 ~ Heat_flow_from_the_earths_core - 1.609,
            0 ~ Heat_gained___needed_for_the_desired_freezing___unfreezing_of_permafrost_ZJ_yr - 0.35770833333333335CH4_in_permafrost_area_melted___frozen_before_heat_constraint_GtC_yr,
            0 ~ Heat_used_in_melting__pos__or_freezing__neg__antarctic_ice_ZJ_yr - 0.0003327Antarctic_ice_losing__pos__or_gaining__neg__GtIce_yr,
            0 ~ Heat_used_in_melting__pos__or_freezing__neg__arctic_sea_ice_ZJ_yr - 8.3175e-7Arctic_ice_melting__pos__or_freezing__neg__km2_yr,
          ]
        end
        function generateEquations11()
          println("#Equation generated:" * "50" * "in: " * "generateEquations11")
          [
            0 ~ Heat_used_in_melting__pos__or_freezing__neg__glacial_ice_ZJ_yr - 0.0003327Glacial_ice_melting__pos__or_freezing__neg__km3_yr,
            0 ~ Heat_used_in_melting__pos__or_freezing__neg__Greenland_ice_that_slid_into_the_water_ZJ_yr - 0.0003327Greenland_ice_that_slid_into_the_ocean_melting__pos__or_freezing__neg__km3_yr,
            0 ~ Heat_used_in_melting__pos__or_freezing__neg__Greenland_ice_ZJ_yr - 0.0003327Greenland_ice_melting__pos__or_freezing__neg__km3_yr,
            0 ~
              Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_antarctic_ice_W_m2 +
              -Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_antarctic_ice_ZJ_yr / UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~ Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_antarctic_ice_ZJ_yr - 0.4Heat_used_in_melting__pos__or_freezing__neg__antarctic_ice_ZJ_yr,
            0 ~
              Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_arctic_ice_W_m2 +
              -Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_arctic_ice_ZJ_yr / UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~ Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_arctic_ice_ZJ_yr - 0.5Heat_used_in_melting__pos__or_freezing__neg__arctic_sea_ice_ZJ_yr,
            0 ~
              Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_Greenland_ice_that_slid_into_the_ocean_W_m2 +
              -Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_Greenland_ice_that_slid_into_the_ocean_ZJ_yr / UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~
              Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_Greenland_ice_that_slid_into_the_ocean_ZJ_yr -
              0.9Heat_used_in_melting__pos__or_freezing__neg__Greenland_ice_that_slid_into_the_water_ZJ_yr,
            0 ~
              Heat_withdrawn_from_atm_by_melting__pos____added__neg__by_freezing_ice_W_m2 +
              -Heat_withdrawn_from_atm_by_melting__pos____added__neg__by_freezing_ice_ZJ_yr / UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~
              (
                (
                  (
                    (
                      (Heat_withdrawn_from_atm_by_melting__pos____added__neg__by_freezing_ice_ZJ_yr - 0.1Heat_used_in_melting__pos__or_freezing__neg__Greenland_ice_that_slid_into_the_water_ZJ_yr) -
                      0.5Heat_used_in_melting__pos__or_freezing__neg__arctic_sea_ice_ZJ_yr
                    ) - 0.6Heat_used_in_melting__pos__or_freezing__neg__antarctic_ice_ZJ_yr
                  ) - Heat_actually_gained___needed_for_freezing___unfreezing_of_permafrost_ZJ_yr
                ) - Heat_used_in_melting__pos__or_freezing__neg__Greenland_ice_ZJ_yr
              ) - Heat_used_in_melting__pos__or_freezing__neg__glacial_ice_ZJ_yr,
            0 ~ (HI_clouds_net_effect__pos_warming__neg_cooling__W_m2 + SW_HI_cloud_efffect_aka_TOA_albedo_W_m2) - LW_HI_cloud_radiation_W_m2,
            0 ~ Hist_Frac_atm_absorption - 0.22058823529411764,
            0 ~
              Human_activity_CH4_emissions -
              ESCIMO_IF_THEN_ELSE(true, CH4_emissions_CO2e_after_exp_12a, ESCIMO_IF_THEN_ELSE(Time > 2020.0, 0.0, CH4_emissions_CO2e_after_exp_12a)),
            0 ~ Human_activity_CH4_emissions_GtCO2e_yr + -Human_activity_CH4_emissions / UNIT_conversion_for_CH4_from_CO2e_to_C,
            0 ~ Humidity_of_atmosphere_current_g_kg - Humidity_of_atmosphere_g_kg,
            0 ~ Humidity_of_atmosphere_g_kg - 0.00125 * Evaporation_as_f_of_temp_ZJ_yr,
            0 ~ Ice_on_land_area_Mkm2 - 1.0e-6 * (Antarctic_ice_area_km2 + Glacial_ice_area_km2 + Greenland_ice_area_km2),
            0 ~
              ((Incoming_solar_W_m2 + ESCIMO_IF_THEN_ELSE(false, 3.0 * ESCIMO_IF_THEN_ELSE((Time > 2015) & (Time < 3.0e7), 1.0, 0.0), 0.0)) - 340.0) -
              Solar_cycle_W_m2,
            0 ~ Incoming_solar_ZJ_yr - Incoming_solar_W_m2 * UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~ InputEmissions_for_tipping_point_search - All_Human_activity_emissions_GtCO2e_yr_Base_for_tipping_point_search,
            0 ~ Intercept_blocked_by_H20_future_equ - 0.33169,
            0 ~ Kyoto_Flour_concentration_ppt - 0.04Kyoto_Flour_gases_in_atm,
            0 ~ Kyoto_Flour_degradation + -Kyoto_Flour_gases_in_atm / Time_to_degrade_Kyoto_Flour_yr,
            0 ~ Kyoto_Flour_emissions - Kyoto_Flour_emissions_after_exp_12a,
            0 ~
              Kyoto_Flour_emissions_after_exp - ESCIMO_IF_THEN_ELSE(
                true,
                Kyoto_Flour_emissions_before_exp,
                ESCIMO_IF_THEN_ELSE(
                  false,
                  Kyoto_Flour_emissions_from_CO2e_C_Roads,
                  ESCIMO_IF_THEN_ELSE(
                    false,
                    Kyoto_Flour_emissions_from_CO2e_CAT,
                    ESCIMO_IF_THEN_ELSE(
                      false,
                      (1.4285714285714286Kyoto_Flour_emissions_pct_contribution_to_Total_CO2e) * Total_CO2e_emissions_as_f_peak__GtCO2e_yr,
                      Kyoto_Flour_emissions_before_exp,
                    ),
                  ),
                ),
              ),
            0 ~
              Kyoto_Flour_emissions_after_exp_12a -
              ESCIMO_IF_THEN_ELSE(true, Kyoto_Flour_emissions_after_exp, ESCIMO_IF_THEN_ELSE(Time > 2020.0, 0.0, Kyoto_Flour_emissions_after_exp)),
            0 ~
              Kyoto_Flour_emissions_before_exp - ESCIMO_IF_THEN_ELSE(
                Time <= 2010,
                Kyoto_Flour_emissions_RCPs_or_JR52,
                ESCIMO_IF_THEN_ELSE(
                  true,
                  Kyoto_Flour_emissions_RCPs_or_JR52,
                  ESCIMO_IF_THEN_ELSE(false, 0.0, ESCIMO_IF_THEN_ELSE(false, Kyoto_Flour_emissions_RCPs_JR_in_2010, 0.0)),
                ),
              ),
            0 ~ Kyoto_Flour_emissions_GtCO2e_yr - 0.006999999999999999Kyoto_Flour_emissions,
            0 ~
              Kyoto_Flour_emissions_RCPs_or_JR52 - ESCIMO_IF_THEN_ELSE(
                false,
                Emissions_of_Kyoto_Flour_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_kt_yr,
                ESCIMO_IF_THEN_ELSE(
                  true,
                  Emissions_of_Kyoto_Flour_with_JR_2052_shape_kt_yr,
                  ESCIMO_IF_THEN_ELSE(
                    false,
                    Emissions_of_Kyoto_Flour_with_JR_2052_shape_kt_yr,
                    ESCIMO_IF_THEN_ELSE(
                      false,
                      Emissions_of_Kyoto_Flour_with_JR_2052_shape_kt_yr,
                      ESCIMO_IF_THEN_ELSE(
                        false,
                        OGHG_Kyoto_Flour_emi_rcp3,
                        ESCIMO_IF_THEN_ELSE(
                          false,
                          OGHG_Kyoto_Flour_emi_rcp45,
                          ESCIMO_IF_THEN_ELSE(false, OGHG_Kyoto_Flour_emi_rcp6, ESCIMO_IF_THEN_ELSE(false, OGHG_Kyoto_Flour_emi_rcp85, 0.0)),
                        ),
                      ),
                    ),
                  ),
                ),
              ),
            0 ~ Land_area_km2 - 1.53e8,
            0 ~ ((Land_covered_with_ice_km2 - Antarctic_ice_area_km2) - Glacial_ice_area_km2) - Greenland_ice_area_km2,
            0 ~ Land_covered_with_ice_Mkm2 - 1.0e-6Land_covered_with_ice_km2,
            0 ~ (LO_clouds_net_effect__pos_warming__neg_cooling__W_m2 + SW_LO_cloud_efffect_aka_cloud_albedo_W_m2) - LW_LO_cloud_radiation_W_m2,
            0 ~ LW_Blocking_multiplier_from_other_GHG - 0.3333333333333333 * (Blocking_multiplier_from_Kyoto_Flour + Blocking_multiplier_from_Montreal_gases + Blocking_multiplier_from_N2O),
            0 ~ LW_Clear_sky_emissions_from_atm - (1.0 - Frac_blocked_by_ALL_GHG) * (BB_radiation_at_Temp_in_atm_ZJ_yr + LW_surface_emissions_escaping_through_atm_window),
            0 ~ LW_Clear_sky_emissions_from_atm_W_m2 + -LW_Clear_sky_emissions_from_atm / UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~ LW_clear_sky_emissions_to_surface - BB_radiation_at_Temp_in_atm_ZJ_yr,
            0 ~ LW_clear_sky_emissions_to_surface_W_m2 + -LW_clear_sky_emissions_to_surface / UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~ (Blocking_of_LW_rad_by_clouds + LW_Cloudy_sky_emissions_from_atm) - LW_Clear_sky_emissions_from_atm,
            0 ~ LW_Cloudy_sky_emissions_from_atm_W_m2 + -LW_Cloudy_sky_emissions_from_atm / UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~ LW_HI_cloud_radiation - LW_HI_cloud_radiation_reference_in_1850_W_m2 * Ratio_of_area_covered_by_high_clouds_current_to_1850 * UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~ LW_HI_cloud_radiation_reference_in_1850_W_m2 - 7.899999999999999,
            0 ~ LW_HI_cloud_radiation_W_m2 + -LW_HI_cloud_radiation / UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~ LW_LO_cloud_radiation - (20.0Ratio_of_area_covered_by_low_clouds_current_to_1850) * UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~ LW_LO_cloud_radiation_W_m2 + -LW_LO_cloud_radiation / UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~ LW_radiation_blocked_by_CH4__pct_ - 100.0Blocked_by_CH4,
            0 ~ LW_radiation_blocked_by_CO2__pct_ - 100.0Blocked_by_CO2,
            0 ~ LW_radiation_blocked_by_H2O__pct_ - 100.0Blocked_by_H20,
            0 ~ LW_radiation_blocked_by_other_GHG__pct_ - 100.0Blocked_by_otherGHG,
          ]
        end
        function generateEquations12()
          println("#Equation generated:" * "50" * "in: " * "generateEquations12")
          [
            0 ~ (LW_re_radiated_by_clouds - LW_HI_cloud_radiation) - LW_LO_cloud_radiation,
            0 ~ LW_re_radiated_by_clouds_W_m2 + -LW_re_radiated_by_clouds / UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~ LW_surface_emission - BB_radiation_at_surface_temp_ZJ_yr,
            0 ~ LW_surface_emission_W_m2 + -LW_surface_emission / UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~ LW_surface_emissions_escaping_through_atm_window - 0.051LW_surface_emission,
            0 ~ LW_surface_emissions_NOT_escaping_through_atm_window - 0.949LW_surface_emission,
            0 ~ LW_surface_emissions_NOT_escaping_through_atm_window_W_m2 + -LW_surface_emissions_NOT_escaping_through_atm_window / UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~ LW_TOA_radiation_from_atm_to_space - LW_Cloudy_sky_emissions_from_atm,
            0 ~ (LW_TOA_radiation_from_atm_to_space + LW_TOA_radiation_from_atm_to_space_difference_wrt_1850) - LW_TOA_radiation_from_atm_to_space_in_1850,
            0 ~ LW_TOA_radiation_from_atm_to_space_difference_wrt_1850_W_m2 + -LW_TOA_radiation_from_atm_to_space_difference_wrt_1850 / UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~ LW_TOA_radiation_from_atm_to_space_W_m2 + -LW_TOA_radiation_from_atm_to_space / UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~ M_2010 - 1720.81,
            0 ~ M_cur - CH4_concentration_ppb,
            0 ~ Man_made_CH4_emissions_pct - ESCIMO_ZIDZ(Human_activity_CH4_emissions, CH4_all_emissions_GtC_yr),
            0 ~ Man_made_fossil_C_emissions_for_cumulation_GtC_yr - Man_made_fossil_C_emissions_GtC_yr,
            0 ~
              Man_made_fossil_C_emissions_GtC_yr -
              ESCIMO_IF_THEN_ELSE(true, CO2_emissions_CO2e_after_exp, ESCIMO_IF_THEN_ELSE(Time > 2020.0, 0.0, CO2_emissions_CO2e_after_exp)),
            0 ~ Man_made_fossil_C_emissions_GtCO2e_yr + -Man_made_fossil_C_emissions_GtC_yr / UNIT_conversion_for_CO2_from_CO2e_to_C,
            0 ~ Melting_constraint_from_the_heat_in__ocean__surface_reservoir - var"combi_Melting_constraint_from_the_heat_in__ocean__surface_reservoir_y[1]",
            0 ~ Melting_constraint_from_the_heat_in_atmosphere_reservoir_fraction - var"combi_Melting_constraint_from_the_heat_in_atmosphere_reservoir_fraction_y[1]",
            0 ~ Melting_restraint_for_permafrost_from_heat_in_atmophere - var"combi_Melting_restraint_for_permafrost_from_heat_in_atmophere_y[1]",
            0 ~ Methane_hydrates_released_and_converted_to_CO2_by_bacteria_GtC - 0.9 * Experimental_release_of_methane,
            0 ~ Methanehydrate_experimental_release_GtC__yr - 0.09999999999999998 * Experimental_release_of_methane,
            0 ~ MODEL_CH4_in_atm_in_ppb - 468.0C_in_atmosphere_in_form_of_CH4,
            0 ~ MODEL_CO2_concentration_in_atmosphere2_ppm + -C_in_atmosphere_GtC / Conversion_constant_GtC_to_ppm,
            0 ~ Model_Volcanic_aerosol_forcing_W_m2 + Volcanic_aerosols_in_stratosphere,
            0 ~ Montreal_emissions_GtCO2e_yr - 0.01Montreal_gases_emissions,
            0 ~ Montreal_gases_concentration_ppt - 0.04Montreal_gases_in_atm,
            0 ~ Montreal_gases_degradation - 0.03333333333333333Montreal_gases_in_atm,
            0 ~ Montreal_gases_emissions - Montreal_gases_emissions_after_exp_12a,
            0 ~
              Montreal_gases_emissions_after_exp_12a - ESCIMO_IF_THEN_ELSE(
                true,
                Montreal_gases_emissions_CO2e_after_exp,
                ESCIMO_IF_THEN_ELSE(Time > 2020.0, 0.0, Montreal_gases_emissions_CO2e_after_exp),
              ),
            0 ~
              Montreal_gases_emissions_before_exp - ESCIMO_IF_THEN_ELSE(
                Time <= 2010,
                Montreal_gases_emissions_RCPs_or_JR52,
                ESCIMO_IF_THEN_ELSE(
                  true,
                  Montreal_gases_emissions_RCPs_or_JR52,
                  ESCIMO_IF_THEN_ELSE(false, 0.0, ESCIMO_IF_THEN_ELSE(false, Montreal_gases_emissions_RCPs_JR_in_2010, 0.0)),
                ),
              ),
            0 ~
              Montreal_gases_emissions_CO2e_after_exp - ESCIMO_IF_THEN_ELSE(
                true,
                Montreal_gases_emissions_before_exp,
                ESCIMO_IF_THEN_ELSE(
                  false,
                  Montreal_gases_emissions_from_CO2e_C_Roads,
                  ESCIMO_IF_THEN_ELSE(
                    false,
                    Montreal_gases_emissions_from_CO2e_CAT,
                    ESCIMO_IF_THEN_ELSE(
                      false,
                      (1.0000000000000002Montreal_gases_emissions_pct_contribution_to_Total_CO2e) * Total_CO2e_emissions_as_f_peak__GtCO2e_yr,
                      Montreal_gases_emissions_before_exp,
                    ),
                  ),
                ),
              ),
            0 ~
              Montreal_gases_emissions_RCPs_or_JR52 - ESCIMO_IF_THEN_ELSE(
                false,
                Emissions_of_Montreal_gases_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_kt_yr,
                ESCIMO_IF_THEN_ELSE(
                  true,
                  Emissions_of_Montreal_gases_with_JR_2052_shape_kt_yr,
                  ESCIMO_IF_THEN_ELSE(
                    false,
                    Emissions_of_Montreal_gases_with_JR_2052_shape_kt_yr,
                    ESCIMO_IF_THEN_ELSE(
                      false,
                      Emissions_of_Montreal_gases_with_JR_2052_shape_kt_yr,
                      ESCIMO_IF_THEN_ELSE(
                        false,
                        OGHG_Montreal_gases_emi_rcp3,
                        ESCIMO_IF_THEN_ELSE(
                          false,
                          OGHG_Montreal_gases_emi_rcp45,
                          ESCIMO_IF_THEN_ELSE(false, OGHG_Montreal_gases_emi_rcp6, ESCIMO_IF_THEN_ELSE(false, OGHG_Montreal_gases_emi_rcp85, 0.0)),
                        ),
                      ),
                    ),
                  ),
                ),
              ),
            0 ~ N_2010 - 363.504,
            0 ~ N_cur - N2O_concentration_ppb,
            0 ~
              N20_emissions_RCPs_or_JR52 - ESCIMO_IF_THEN_ELSE(
                false,
                Emissions_of_anthro_N20_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_MtN2O_yr,
                ESCIMO_IF_THEN_ELSE(
                  true,
                  Emissions_of_anthro_N20_with_JR_2052_shape_MtN20_yr,
                  ESCIMO_IF_THEN_ELSE(
                    false,
                    Emissions_of_anthro_N20_with_JR_2052_shape_MtN20_yr,
                    ESCIMO_IF_THEN_ELSE(
                      false,
                      Emissions_of_anthro_N20_with_JR_2052_shape_MtN20_yr,
                      ESCIMO_IF_THEN_ELSE(
                        false,
                        othGHG_N20_man_made_emissions_rcp3,
                        ESCIMO_IF_THEN_ELSE(
                          false,
                          othGHG_N20_man_made_emissions_rcp45,
                          ESCIMO_IF_THEN_ELSE(
                            false,
                            othGHG_N20_man_made_emissions_rcp6,
                            ESCIMO_IF_THEN_ELSE(false, othGHG_N20_man_made_emissions_rcp85, 0.0),
                          ),
                        ),
                      ),
                    ),
                  ),
                ),
              ),
            0 ~ N2O_concentration_ppb - 0.305N2O_in_atmosphere_MtN2O,
            0 ~ N2O_degradation_MtN2O_yr - 0.010526315789473684N2O_in_atmosphere_MtN2O,
            0 ~
              N2O_man_made_emissions - ESCIMO_IF_THEN_ELSE(
                Time <= 2010,
                N20_emissions_RCPs_or_JR52,
                ESCIMO_IF_THEN_ELSE(
                  true,
                  N20_emissions_RCPs_or_JR52,
                  ESCIMO_IF_THEN_ELSE(false, 0.0, ESCIMO_IF_THEN_ELSE(false, N20_emissions_RCPs_JR_in_2010, 0.0)),
                ),
              ),
            0 ~
              N2O_man_made_emissions_after_exp - ESCIMO_IF_THEN_ELSE(
                true,
                N2O_man_made_emissions,
                ESCIMO_IF_THEN_ELSE(
                  false,
                  N2O_man_made_emissions_from_CO2e_C_Roads,
                  ESCIMO_IF_THEN_ELSE(
                    false,
                    N2O_man_made_emissions_from_CO2e_CAT,
                    ESCIMO_IF_THEN_ELSE(
                      false,
                      (0.03355704697986577N2O_emissions_pct_contribution_to_Total_CO2e) * Total_CO2e_emissions_as_f_peak__GtCO2e_yr,
                      N2O_man_made_emissions,
                    ),
                  ),
                ),
              ),
            0 ~
              N2O_man_made_emissions_exp_12a -
              ESCIMO_IF_THEN_ELSE(true, N2O_man_made_emissions_after_exp, ESCIMO_IF_THEN_ELSE(Time > 2020.0, 0.0, N2O_man_made_emissions_after_exp)),
            0 ~ N2O_man_made_emissions_GtCO2e_yr - 0.298N2O_man_made_emissions_exp_12a,
            0 ~ NatEvent_d__slowing_down_ocean_circulation_from_2015 - ESCIMO_IF_THEN_ELSE(Time > 2015, 1.0, 1.0),
            0 ~ (Natural_CH4_emissions - 0.19) - CH4_emissions_from_wetlands_destruction,
            0 ~ Natural_CH4_emissions_pct - ESCIMO_ZIDZ(Natural_CH4_emissions, CH4_all_emissions_GtC_yr),
            0 ~ NATURE_CCS_Fig3_GtC_yr,
            0 ~ NATURE_CCS_removal_experiment_multiplier - var"combi_NATURE_CCS_removal_experiment_multiplier_y[1]",
            0 ~ (2.0 + Net_additions_to_C_in_TUNDRA_DeadB_and_soil_GtC) - C_in_TUNDRA_DeadB_and_soil_GtC,
            0 ~ (2.0 + Net_additions_to_C_in_TUNDRA_LB_GtC) - C_in_TUNDRA_LB_GtC,
            0 ~ (Flow_from_biomass_to_atm_Gtc_pr_yr + Net_C_flow_from_atm_to_biomass_GtC_pr_yr) - Flow_from_atm_to_biomass_GtC_pr_yr,
          ]
        end
        function generateEquations13()
          println("#Equation generated:" * "50" * "in: " * "generateEquations13")
          [
            0 ~
              (
                (
                  (
                    (
                      (
                        (
                          (
                            CO2_flux_from_atm_to_GRASS_for_new_growth_GtC_yr +
                            CO2_flux_from_atm_to_NF_for_new_growth_GtC_yr +
                            CO2_flux_from_atm_to_TROP_for_new_growth_GtC_yr +
                            CO2_flux_from_atm_to_TUNDRA_for_new_growth_GtC_yr +
                            C_diffusion_into_ocean_from_atm +
                            Carbon_captured_and_stored_GtC___yr +
                            Net_C_to_atm
                          ) - Avg_volcanic_activity_GtC_yr
                        ) - CH4_in_the_atmosphere_converted_to_CO2
                      ) - CO2_flux_GRASS_to_atm_Gtc_yr
                    ) - CO2_flux_NF_to_atm_Gtc_yr
                  ) - CO2_flux_TROP_to_atm_GtC_yr
                ) - CO2_flux_TUNDRA_to_atm_Gtc_yr
              ) - Man_made_fossil_C_emissions_GtC_yr,
            0 ~ Net_C_to_atm_rate - Net_C_to_atm,
            0 ~ (CO2_flux_GRASS_to_atm_Gtc_yr + Net_CO2_flow_between_grass_and_atmosphere_GtC) - CO2_flux_from_atm_to_GRASS_for_new_growth_GtC_yr,
            0 ~ (CO2_flux_TUNDRA_to_atm_Gtc_yr + Net_CO2_flow_between_TUNDRA_and_atmosphere_GtC) - CO2_flux_from_atm_to_TUNDRA_for_new_growth_GtC_yr,
            0 ~
              (
                (
                  (
                    Convection_aka_sensible_heat_flow +
                    Evaporation_aka_latent_heat_flow +
                    Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_Greenland_ice_that_slid_into_the_ocean_ZJ_yr +
                    Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_antarctic_ice_ZJ_yr +
                    Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_arctic_ice_ZJ_yr +
                    LW_surface_emission +
                    Net_flow_of_heat_into_surface +
                    Net_heat_flow_ocean_from_surface_to_deep__ZJ_yr_
                  ) - LW_clear_sky_emissions_to_surface
                ) - LW_re_radiated_by_clouds
              ) - SW_surface_absorption,
            0 ~ ((Biological_removal_of_C_from_WSW_GtC_per_yr + Depositing_of_C_to_sediment + Net_flux_to_ocean_GtC_yr) - C_diffusion_into_ocean_from_atm) - C_runoff_from_biomass_soil,
            0 ~ Net_heat_flow_ocean_between_surface_and_deep_per_K_of_difference_ZJ_yr_K - ESCIMO_IF_THEN_ELSE(Time > 2020, 10.0, 10.0),
            0 ~ Net_heat_flow_ocean_from_surface_to_deep__ZJ_yr_ - Net_heat_flow_ocean_between_surface_and_deep_per_K_of_difference_ZJ_yr_K * Surface_deep__ocean__temp_diff_degC,
            0 ~ Net_heat_flow_ocean_from_surface_to_deep_W_m2 + -Net_heat_flow_ocean_from_surface_to_deep__ZJ_yr_ / UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~
              (
                (
                  (
                    (
                      Heat_withdrawn_from_atm_by_melting__pos____added__neg__by_freezing_ice_ZJ_yr +
                      LW_TOA_radiation_from_atm_to_space +
                      LW_clear_sky_emissions_to_surface +
                      Net_heat_flow_to_atm_ZJ_yr__needed_for_comparisons_with_history_
                    ) - Convection_aka_sensible_heat_flow
                  ) - Evaporation_aka_latent_heat_flow
                ) - LW_surface_emissions_NOT_escaping_through_atm_window
              ) - SW_Atmospheric_absorption,
            0 ~ Net_marine_primary_production_NMPP_GtC_pr_yr - (0.4 * Effect_of_C_concentration_on_NMPP) * Effect_of_acidification_on_NMPP * Effect_of_temperature_on_NMPP,
            0 ~ (NEW_Temp_ocean_surface_in_1850_in_K - 273.15) - Temp__ocean__surface_in_1850_C,
            0 ~ NF_Avg_life_biomass_yr - ESCIMO_IF_THEN_ELSE(Time > 2020, 60.0, 60.0),
            0 ~ NF_being_deforested_Mkm2_yr - NF_historical_deforestation_pct_yr * NF_with_normal_cover_Mkm2,
            0 ~ NF_being_harvested_by_clear_cutting_Mkm2_yr - NF_being_harvested_Mkm2_yr * NF_clear_cut_fraction,
            0 ~
              NF_being_harvested_Mkm2_yr +
              ((-1000.0NF_usage_cutoff) * Use_of_NF_biomass_for_energy_GtBiomass_yr * ESCIMO_IF_THEN_ELSE(false, POLICY_4_Stopping_logging_in_Northern_forests, 1.0)) /
              NF_living_biomass_densitiy_tBiomass_pr_km2,
            0 ~ NF_being_harvested_normally_Mkm2_yr - NF_being_harvested_Mkm2_yr * (1.0 - NF_clear_cut_fraction),
            0 ~
              NF_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr -
              (0.001NF_living_biomass_densitiy_tBiomass_pr_km2) *
              (NF_being_deforested_Mkm2_yr + NF_being_harvested_by_clear_cutting_Mkm2_yr + NF_being_harvested_normally_Mkm2_yr + NF_burning_Mkm2_yr),
            0 ~ NF_Biomass_in_construction_material_being_burnt_GtBiomass_yr - 0.025NF_Biomass_locked_in_construction_material_GtBiomass,
            0 ~ NF_biomass_new_growing_GtBiomass___yr + -NF_potential_less_actual_living_biomass_GtBiomass / NF_Speed_of_regrowth_yr,
            0 ~ NF_burning_Mkm2_yr - (0.006999999999999999 * Effect_of_temperature_on_fire_incidence_dimensionless) * NF_with_normal_cover_Mkm2,
            0 ~ NF_clear_cut_fraction - var"combi_NF_clear_cut_fraction_y[1]",
            0 ~ NF_Dead_biomass_decomposing_GtBiomass_yr - 0.004NF_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass,
            0 ~ NF_DeadB_and_SOM_tB_per_km2 - (27500.0 * Effect_of_CO2_on_new_biomass_growth) * Effect_of_temperature_on_new_biomass_growth_dimensionless,
            0 ~ NF_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr - (0.001NF_DeadB_and_SOM_tB_per_km2) * NF_being_deforested_Mkm2_yr,
            0 ~ NF_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr - (0.0001NF_DeadB_and_SOM_tB_per_km2) * NF_being_harvested_normally_Mkm2_yr,
            0 ~ NF_for_construction_use_GtBiomass_yr - Use_of_NF_biomass_for_construction_GtBiomass_yr,
            0 ~ NF_historical_deforestation_pct_yr - 0.0002,
            0 ~ NF_land_taken_out_of_use_GtBiomass - (0.001NF_land_taken_out_of_use_Mkm2) * NF_living_biomass_densitiy_tBiomass_pr_km2,
            0 ~ (((NF_land_taken_out_of_use_Mkm2 - NF_area_burnt_Mkm2) - NF_area_clear_cut_Mkm2) - NF_area_deforested_Mkm2) - NF_area_harvested_Mkm2,
            0 ~ NF_living_biomass_densitiy_tBiomass_pr_km2 - (7500.0 * Effect_of_CO2_on_new_biomass_growth) * Effect_of_temperature_on_new_biomass_growth_dimensionless,
            0 ~ NF_Living_biomass_rotting_GtBiomass_yr + -NF_Living_biomass_GtBiomass / NF_Avg_life_biomass_yr,
            0 ~ (NF_Living_biomass_GtBiomass + NF_potential_less_actual_living_biomass_GtBiomass) - NF_potential_living_biomass_GtBiomass,
            0 ~ NF_potential_living_biomass_GtBiomass - (0.001NF_living_biomass_densitiy_tBiomass_pr_km2) * (NF_potential_area_Mkm2 - NF_area_deforested_Mkm2),
            0 ~ NF_regrowing_after_being_burnt_Mkm2_yr + -NF_area_burnt_Mkm2 / Time_to_regrow_NF_after_buning_yr,
            0 ~ NF_regrowing_after_being_clear_cut_Mkm2_yr + -NF_area_clear_cut_Mkm2 / (2.0Time_to_regrow_NF_after_buning_yr),
            0 ~ NF_regrowing_after_being_deforested_Mkm2_yr - 0.0125NF_area_deforested_Mkm2,
            0 ~ NF_regrowing_after_harvesting_Mkm2_yr + -NF_area_harvested_Mkm2 / Time_to_regrow_NF_after_buning_yr,
            0 ~ NF_runoff - 0.0005NF_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass,
            0 ~ NF_soil_degradation_from_clear_cutting_GtBiomass_yr - (0.0005NF_DeadB_and_SOM_tB_per_km2) * NF_being_harvested_by_clear_cutting_Mkm2_yr,
            0 ~ NF_soil_degradation_from_forest_fires_GtBiomass_yr,
            0 ~ NF_Speed_of_regrowth_yr - ESCIMO_IF_THEN_ELSE(Time > 2020, 3.0, 3.0),
            0 ~
              (
                (
                  (
                    (
                      (NF_Sum_outflows_GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass_yr - NF_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr) -
                      NF_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr
                    ) - NF_Dead_biomass_decomposing_GtBiomass_yr
                  ) - NF_runoff
                ) - NF_soil_degradation_from_clear_cutting_GtBiomass_yr
              ) - NF_soil_degradation_from_forest_fires_GtBiomass_yr,
            0 ~ NF_TUNDRA_Biomass_in_construction_material_left_to_rot_GtBiomass_yr - 0.025NF_Biomass_locked_in_construction_material_GtBiomass,
            0 ~ NF_usage_as_pct_of_potial_area + (((-NF_area_burnt_Mkm2 - NF_area_clear_cut_Mkm2) - NF_area_deforested_Mkm2) - NF_area_harvested_Mkm2) / NF_with_normal_cover_Mkm2,
            0 ~ NF_usage_cutoff - var"combi_NF_usage_cutoff_y[1]",
            0 ~ (NF_area_burnt_Mkm2 + NF_area_clear_cut_Mkm2 + NF_area_deforested_Mkm2 + NF_area_harvested_Mkm2 + NF_with_normal_cover_Mkm2) - NF_potential_area_Mkm2,
            0 ~ Ocean_area_km2 - 3.57e8,
            0 ~ Ocean_circulation_slowdown_from_Greenland_ice_sliding_into_the_Atlantic - max(0.6699999999999999, min(1.0, 1.0 - 0.007444444444444444Time_less_Greenland_slide_experiment_start_yr)),
            0 ~
              Ocean_heat_used_for_melting_ZJ_yr +
              (
                (
                  -Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_Greenland_ice_that_slid_into_the_ocean_ZJ_yr -
                  Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_antarctic_ice_ZJ_yr
                ) - Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_arctic_ice_ZJ_yr
              ) / Heat_in_surface,
          ]
        end
        function generateEquations14()
          println("#Equation generated:" * "50" * "in: " * "generateEquations14")
          [
            0 ~ Ocean_surface_area_km2 - 3.619e8,
            0 ~ (Ocean_surface_delta_temp_to_1850_C + Temp__ocean__surface_in_1850_C) - Temp__ocean__surface_in_K,
            0 ~ (Open_water_as_frac_of_ocean_area + Arctic_ice__on_sea__area_km2 / Ocean_area_km2) - 1.0,
            0 ~ ((LW_HI_cloud_radiation_W_m2 + Outgoing_radiation_at_TOA_W_m2) - LW_TOA_radiation_from_atm_to_space_W_m2) - Reflected_Solar_SW_W_m2,
            0 ~ pct_change_in_fraction_blocked_by_ALL_GHG_wrt_1850 + (-100.0 * (Frac_blocked_by_ALL_GHG - Fraction_blocked_by_ALL_GHG_in_1850)) / Fraction_blocked_by_ALL_GHG_in_1850,
            0 ~ pct_change_in_fraction_blocked_by_C02_wrt_1850 + (-100.0 * (Blocked_by_CO2 - Fraction_blocked_CO2_in_1850)) / Fraction_blocked_CO2_in_1850,
            0 ~ pct_change_in_fraction_blocked_by_CH4_wrt_1850 + (-100.0 * (Blocked_by_CH4 - Fraction_blocked_CH4_in_1850)) / Fraction_blocked_CH4_in_1850,
            0 ~ pct_change_in_fraction_blocked_by_othGHG_wrt_1850 + (-100.0 * (Blocked_by_otherGHG - Fraction_blocked_othGHG_in_1850)) / Fraction_blocked_othGHG_in_1850,
            0 ~ pct_reduction_in_C_in_GRASS + (-100.0 * (init_C_in_GRASS - C_in_GRASS_GtC)) / init_C_in_GRASS,
            0 ~ pct_reduction_in_C_in_NF + (-100.0 * (init_C_in_NF - C_in_NF_GtC)) / init_C_in_NF,
            0 ~ pct_reduction_in_C_in_TROP + (-100.0 * (init_C_in_TROP - C_in_TROP_GtC)) / init_C_in_TROP,
            0 ~ pct_reduction_in_C_in_TUNDRA + (-100.0 * (init_C_in_TUNDRA - C_in_TUNDRA_GtC)) / init_C_in_TUNDRA,
            0 ~ Permafrost_area_km2 - 20833.333333333332C_in_permafrost_in_form_of_CH4,
            0 ~ Permafrost_CH4_emissions_pct - ESCIMO_ZIDZ(CH4_release___capture_from_permafrost_area_loss___gain_GtC_yr, CH4_all_emissions_GtC_yr),
            0 ~ Permafrost_melting_cutoff - var"combi_Permafrost_melting_cutoff_y[1]",
            0 ~ pH_in_cold_deep_water + (-163.2 * (0.9997 - 0.0017Temp__ocean__deep_in_C)) / CC_in_deep_box_ymoles_per_litre__dimensionless_^0.385,
            0 ~ ph_in_cold_downwelling_water + (-163.2 * (0.9997 - 0.0017Temp_of_cold_downwelling_water)) / CC_in_cold_downwelling_ymoles_per_litre__dimensionless_^0.385,
            0 ~ pH_in_cold_suface_water + (-163.2 * (0.9997 - 0.0017Temp_of_cold_surface_water)) / CC_in_cold_surface_ymoles_per_litre__dimensionless_^0.385,
            0 ~
              pH_in_surface +
              (-Volume_cold_ocean_0_to_100m * pH_in_cold_suface_water) / (Volume_cold_ocean_0_to_100m + Volume_warm_ocean_0_to_100m) +
              (-Volume_warm_ocean_0_to_100m * pH_in_warm_surface_water) / (Volume_cold_ocean_0_to_100m + Volume_warm_ocean_0_to_100m),
            0 ~ pH_in_upwelling_water + (-163.2 * (0.9997 - 0.00085 * (Temp__ocean__deep_in_C + Temp_surface_C))) / CC_in_intermediate_box_ymoles_per_litre__dimensionless_^0.385,
            0 ~ pH_in_warm_surface_water + (-163.2 * (0.9997 - 0.0017Temp_surface_C)) / CC_in_warm_surface_ymoles_per_litre__dimensionless_^0.385,
            0 ~ POLICY_4_Stopping_logging_in_Northern_forests - ESCIMO_IF_THEN_ELSE(Time > 2015, 0.0, 1.0),
            0 ~ (Outgoing_radiation_at_TOA_W_m2 + Radiation_balance_at_TOA_in_less_out_W_m2) - Incoming_solar_W_m2,
            0 ~ Radiative_forcing_from_CH4_wrt_1850_W_m2 - Blocked_by_CH4 * RF_from_Ga__BB_radiation_less_TOA_radiation_W_m2,
            0 ~ Radiative_forcing_from_CO2_wrt_1850_W_m2 - Blocked_by_CO2 * RF_from_Ga__BB_radiation_less_TOA_radiation_W_m2,
            0 ~ Radiative_forcing_from_H2O_wrt_1850_W_m2 - Blocked_by_H20 * RF_from_Ga__BB_radiation_less_TOA_radiation_W_m2,
            0 ~ Radiative_forcing_from_othGHG_wrt_1850_W_m2 - Blocked_by_otherGHG * RF_from_Ga__BB_radiation_less_TOA_radiation_W_m2,
            0 ~ (LW_Clear_sky_emissions_from_atm_W_m2 + Radiative_forcing_wrt_1850_W_m2_0) - 2.0,
            0 ~ (Rate_of_destruction_of_wetlands + ESCIMO_STEP(t, 0.0, 2025.0)) - ESCIMO_STEP(t, 0.0, 2020.0),
            0 ~ Ratio_of_area_covered_by_high_clouds_current_to_1850 - 5.0Area_covered_by_high_clouds,
            0 ~ Ratio_of_area_covered_by_low_clouds_current_to_1850 - 2.5Area_covered_by_low_clouds,
            0 ~ RCPFossil_fuel_usage_cutoff - var"combi_RCPFossil_fuel_usage_cutoff_y[1]",
            0 ~ (((Reflected_Solar_SW - SW_HI_cloud_efffect_aka_cloud_albedo) - SW_LO_cloud_efffect_aka_cloud_albedo) - SW_clear_sky_reflection_aka_scattering) - SW_surface_reflection,
            0 ~ Reflected_Solar_SW_W_m2 + -Reflected_Solar_SW / UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~ ((RF_CH4_IPCC_formula_W_m2 + f_M_cur_N_2010_) - 0.0594 * (sqrt(CH4_concentration_ppb) - 1.4142135623730951)) - f_M_1750_N_2010__for_ch4_forcing,
            0 ~ RF_CO2_Model_Myhre_formula - 5.35 * ESCIMO_ln(0.5CO2_concentration_used__after_any_experiments__ppm),
            0 ~ RF_CO2_Model_Myhre_formula_1850 - 5.35 * ESCIMO_ln(CO2_concentration_used__after_any_experiments__ppm / CO2_concentration_in_1850_ppm),
            0 ~ RF_CO2_RCP3_Myhre_formula - 5.35 * ESCIMO_ln(0.5RCP_3_CO2_concentration_1850_2100_ppm),
            0 ~ RF_CO2_RCP45_Myhre_formula - 5.35 * ESCIMO_ln(0.5RCP_45_CO2_concentration_1850_2100_ppm),
            0 ~ RF_CO2_RCP6_Myhre_formula - 5.35 * ESCIMO_ln(0.5RCP_6_CO2_concentration_1850_2100_ppm),
            0 ~ RF_CO2_RCP85_Myhre_formula - 5.35 * ESCIMO_ln(0.5RCP_85_CO2_concentration_1850_2100_ppm),
            0 ~ (2.0 + RF_from_Ga__BB_radiation_less_TOA_radiation_W_m2) - Ga__BB_radiation_less_TOA_radiation_W_m2,
            0 ~ ((RF_N20_IPCC_formula_W_m2 + f_M_2010_N_cur_) - 0.12 * (sqrt(N2O_concentration_ppb) - 1.4142135623730951)) - f_M2010_N_1750__for_n20_forcing,
            0 ~ (Sea_level_change_from_melting_ice_and_thermal_expansion_m - Sea_level_rise_from_melting_ice_m) - Total_sea_level_change_from_thermal_expansion_m,
            0 ~ Sea_level_change_from_thermal_expansion_deep_m + -Volume_expansion_from_thermal_expansion_deep_Gm3_km3 / Ocean_surface_area_km2,
            0 ~ Sea_level_change_from_thermal_expansion_surface_m + -Volume_expansion_from_thermal_expansion_surface_Gm3_km3 / Ocean_surface_area_km2,
            0 ~ Sea_level_rise_from_melting_ice_m + (-1000.0Cumulative_ocean_volume_increase_due_to_ice_melting_km3) / Ocean_surface_area_km2,
            0 ~ Sea_level_rise_history_m - Sea_level_rise_history_mm * UNIT_conversion_mm_to_m,
            0 ~ Seconds_per_yr - 3.1536e7,
            0 ~ Sensitivity_of_high_cloud_coverage_to_temp - ESCIMO_IF_THEN_ELSE(Time > 2020, 50.0, 50.0),
          ]
        end
        function generateEquations15()
          println("#Equation generated:" * "50" * "in: " * "generateEquations15")
          [
            0 ~
              Shifting_GRASS_to_DESERT_Mkm2_yr - ESCIMO_IF_THEN_ELSE(
                Temp_driver_to_shift_biomes_degC > 0.0,
                ((0.1GRASS_potential_area_Mkm2) * Temp_driver_to_shift_biomes_degC^2) / Effect_of_humidity_on_shifting_biomes,
                0.0,
              ),
            0 ~
              Shifting_GRASS_to_NF_Mkm2_yr +
              ESCIMO_IF_THEN_ELSE(Temp_driver_to_shift_biomes_degC < 0.0, (0.002GRASS_potential_area_Mkm2) * Temp_driver_to_shift_biomes_degC, 0.0),
            0 ~
              Shifting_GRASS_to_TROP_Mkm2_yr +
              ESCIMO_IF_THEN_ELSE(Temp_driver_to_shift_biomes_degC < 0.0, (0.004GRASS_potential_area_Mkm2) * Temp_driver_to_shift_biomes_degC, 0.0),
            0 ~ Shifting_ice_on_land_to_tundra_Mkm2_yr - Shifting_ice_to_tundra_from_detail_ice_on_land_Mkm2_pr_yr,
            0 ~ ((Shifting_ice_to_tundra_from_detail_ice_on_land_Mkm2_pr_yr - Antarctic_ice_area_decrease_Mkm2_pr_yr) - Glacial_ice_area_decrease_Mkm2_pr_yr) - Greenland_ice_area_decrease_Mkm2_pr_yr,
            0 ~
              Shifting_NF_to_GRASS_Mkm2_yr -
              ESCIMO_IF_THEN_ELSE(Temp_driver_to_shift_biomes_degC > 0.0, (0.0002NF_potential_area_Mkm2) * Temp_driver_to_shift_biomes_degC, 0.0),
            0 ~
              Shifting_NF_to_TROP_Mkm2_yr - ESCIMO_IF_THEN_ELSE(Temp_driver_to_shift_biomes_degC > 0.0, (0.004NF_potential_area_Mkm2) * Temp_driver_to_shift_biomes_degC, 0.0),
            0 ~
              Shifting_NF_to_Tundra_Mkm2_yr +
              ESCIMO_IF_THEN_ELSE(Temp_driver_to_shift_biomes_degC < 0.0, (0.002NF_potential_area_Mkm2) * Temp_driver_to_shift_biomes_degC, 0.0),
            0 ~
              Shifting_TROP_to_GRASS_Mkm2_yr - ESCIMO_IF_THEN_ELSE(
                Temp_driver_to_shift_biomes_degC > 0.0,
                ((0.001TROP_potential_area_Mkm2) * Temp_driver_to_shift_biomes_degC) / Effect_of_humidity_on_shifting_biomes,
                0.0,
              ),
            0 ~
              Shifting_TROP_to_NF_Mkm2_yr +
              ESCIMO_IF_THEN_ELSE(Temp_driver_to_shift_biomes_degC < 0.0, (0.02TROP_potential_area_Mkm2) * Temp_driver_to_shift_biomes_degC, 0.0),
            0 ~ ((Shifting_tundra_to_ice_from_detail_ice_on_land_Mkm2_pr_yr - Antarctic_ice_area_increase_Mkm2_pr_yr) - Glacial_ice_area_increase_Mkm2_pr_yr) - Greenland_ice_area_increase_Mkm2_pr_yr,
            0 ~ Shifting_tundra_to_ice_on_land_Mkm2_yr - Shifting_tundra_to_ice_from_detail_ice_on_land_Mkm2_pr_yr,
            0 ~
              Shifting_Tundra_to_NF_Mkm2_yr -
              ESCIMO_IF_THEN_ELSE(Temp_driver_to_shift_biomes_degC > 0.0, (0.004Temp_driver_to_shift_biomes_degC) * Tundra_potential_area_Mkm2, 0.0),
            0 ~ SHUT_OFF_permafrost - ESCIMO_IF_THEN_ELSE(Time > 500000.0, 0.0, 1.0),
            0 ~
              Sifting_DESERT_to_GRASS_Mkm2_yr + ESCIMO_IF_THEN_ELSE(
                Temp_driver_to_shift_biomes_degC < 0.0,
                (0.02DESERT_Mkm2) * Slope_of_effect_of_temp_shifting_DESERT_to_GRASS * Temp_driver_to_shift_biomes_degC,
                0.0,
              ),
            0 ~ Slider_for_H2O_slope - ESCIMO_IF_THEN_ELSE(Time > 2020, 0.0, 0.0),
            0 ~ Slope_blocked_by_H20_future_equ - 0.21,
            0 ~ Slope_btw_temp_and_permafrost_melting___freezing - ESCIMO_IF_THEN_ELSE(Time < 2020, 1.0, 1.0),
            0 ~ Slope_of_effect_of_temp_shifting_DESERT_to_GRASS - ESCIMO_IF_THEN_ELSE(Time > 2020, 0.4, 0.4),
            0 ~ Slope_temp_vs_glacial_ice_melting - ESCIMO_IF_THEN_ELSE(Time > 2020, 1.0, 1.0),
            0 ~ Slowing_of_recapture_of_CH4_dmnl - ESCIMO_IF_THEN_ELSE(Effect_of_temp_on_permafrost_melting_dmnl < 0.0, 0.01, 1.0),
            0 ~ Snowball_earth_cutoff - var"combi_Snowball_earth_cutoff_y[1]",
            0 ~ Solar_cycle_W_m2 - ESCIMO_IF_THEN_ELSE(Time > 2011, Solar_sine_forcing_W_m2, Historical_forcing_from_solar_insolation_W_m2),
            0 ~ (Solar_sine_forcing_W_m2 - 0.05) - (0.1 * ESCIMO_IF_THEN_ELSE(Time > 2015, 1.0, 1.0)) * sin(0.5709090909090909 * (3.5 + Time)),
            0 ~ Stop_of_human_deforestation - ESCIMO_IF_THEN_ELSE(Time > 3000.0, 0.0, 1.0),
            0 ~ (((((Sum_biomes_Mkm2 - 1.0e-6Land_covered_with_ice_km2) - DESERT_Mkm2) - GRASS_potential_area_Mkm2) - NF_potential_area_Mkm2) - TROP_potential_area_Mkm2) - Tundra_potential_area_Mkm2,
            0 ~ (((sum_blocked - Blocked_by_CH4) - Blocked_by_CO2) - Blocked_by_H20) - Blocked_by_otherGHG,
            0 ~ (Sum_heat_to_ocean_1972_to_2008_ZJ - Sum_heat_to_deep_ocean_btw_72_and_08) - Sum_heat_to_surface_ocean_btw_72_and_08,
            0 ~
              (
                (
                  (
                    (Sum_outflows_GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass_yr - GRASS_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr) -
                    GRASS_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr
                  ) - GRASS_Dead_biomass_decomposing_GtBiomass_yr
                ) - GRASS_runoff
              ) - GRASS_soil_degradation_from_forest_fires_GtBiomass_yr,
            0 ~ (Surface_deep__ocean__temp_diff_degC + Temp_ocean_deep_in_K) - Temp__ocean__surface_in_K,
            0 ~
              (
                ((Convection_aka_sensible_heat_flow + Evaporation_aka_latent_heat_flow + LW_surface_emission + Surface_imbalance_pos_is_TO_surface) - LW_clear_sky_emissions_to_surface) -
                LW_re_radiated_by_clouds
              ) - SW_surface_absorption,
            0 ~ Surface_imbalance_pos_is_TO_surface_W_m2 + -Surface_imbalance_pos_is_TO_surface / UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~ (Surface_ocean__warm__volume - Intermediate_upwelling_water_volume_100m_to_1km_Gm3) - Warm_surface_water_volume_Gm3,
            0 ~ SW_Atmospheric_absorption - Frac_atm_absorption * Incoming_solar_ZJ_yr,
            0 ~ SW_Atmospheric_absorption_W_m2 + -SW_Atmospheric_absorption / UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~ (SW_clear_sky_reflection_aka_scattering + Total_net_aerosol_forcing_ZJ_yr) - 0.0837Incoming_solar_ZJ_yr,
            0 ~ SW_clear_sky_reflection_aka_scattering_W_m2 + -SW_clear_sky_reflection_aka_scattering / UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~ SW_HI_cloud_efffect_aka_cloud_albedo - (0.006Incoming_solar_ZJ_yr) * Ratio_of_area_covered_by_high_clouds_current_to_1850,
            0 ~ SW_HI_cloud_efffect_aka_TOA_albedo_W_m2 + -SW_HI_cloud_efffect_aka_cloud_albedo / UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~ SW_LO_cloud_efffect_aka_cloud_albedo - (0.158Incoming_solar_ZJ_yr) * Ratio_of_area_covered_by_low_clouds_current_to_1850,
            0 ~ SW_LO_cloud_efffect_aka_cloud_albedo_W_m2 + -SW_LO_cloud_efffect_aka_cloud_albedo / UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~ (SW_surface_absorption + SW_surface_reflection) - SW_to_surface,
            0 ~ (2.0 + SW_surface_absorption_W_m2_wrt_1850) - SW_surface_absorption_W_m2,
            0 ~ SW_surface_absorption_W_m2 + -SW_surface_absorption / UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~ SW_surface_reflection - Avg_earths_surface_albedo * SW_to_surface,
            0 ~ (2.0 + SW_surface_reflection_W_m2_wrt_1850) - SW_surface_reflection_W_m2,
            0 ~ SW_surface_reflection_W_m2 + -SW_surface_reflection / UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~
              (SW_Atmospheric_absorption + SW_HI_cloud_efffect_aka_cloud_albedo + SW_LO_cloud_efffect_aka_cloud_albedo + SW_clear_sky_reflection_aka_scattering + SW_to_surface) - Incoming_solar_ZJ_yr,
            0 ~ SW_to_surface_W_m2 + -SW_to_surface / UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~ Temp__ocean__deep_in_1850_in_K - 277.15,
          ]
        end
        function generateEquations16()
          println("#Equation generated:" * "50" * "in: " * "generateEquations16")
          [
            0 ~ (273.15 + Temp__ocean__deep_in_C) - Temp_ocean_deep_in_K,
            0 ~ (9.7 + Temp__ocean__surface_in_K) - Temp_surface_average_K,
            0 ~ Temp_atm_average_K - Conversion_heat_atm_to_temp * Heat_in_atmosphere_ZJ,
            0 ~ (273.15 + Temp_atm_in_C) - Temp_atm_average_K,
            0 ~ Temp_driver_to_shift_biomes_degC - 0.0002679101448128004Temp_surface_C,
            0 ~ Temp_gradient - 0.25Temp_diff_relevant_for_melting_or_freezing_from_1850,
            0 ~ (1.0 + Temp_gradient_minus_1) - Temp_gradient,
            0 ~ Temp_gradient_minus_1___slope - Slope_btw_temp_and_permafrost_melting___freezing * Temp_gradient_minus_1,
            0 ~ Temp_ocean_deep_in_K - Conversion_constant_heat_ocean_deep_to_temp * Heat_in_deep_ZJ,
            0 ~ Temp_of_cold_downwelling_water - 0.5 * (Temp__ocean__deep_in_C + Temp_of_cold_surface_water),
            0 ~ Temp_of_cold_surface_water - 0.3333333333333333Temp_surface_C,
            0 ~ (13.66500000000002 + Temp_surface_anomaly_compared_to_1850_degC) - Temp_surface_C,
            0 ~ Temp_surface_average_K - Conversion_heat_surface_to_temp * Heat_in_surface,
            0 ~ (273.15 + Temp_surface_C) - Temp_surface_average_K,
            0 ~ Temp_surface_current_divided_by_value_in_1850_K_K - 0.0034865679967923573Temp_surface_average_K,
            0 ~ Thermal_expansion_deep_in_1850_pct - var"combi_Thermal_expansion_deep_in_1850_pct_y[1]",
            0 ~ Thermal_expansion_deep_pct - var"combi_Thermal_expansion_deep_pct_y[1]",
            0 ~ Thermal_expansion_surface_in_1850_pct - var"combi_Thermal_expansion_surface_in_1850_pct_y[1]",
            0 ~ Thermal_expansion_surface_pct - var"combi_Thermal_expansion_surface_pct_y[1]",
            0 ~ Time_in_trunk - ESCIMO_IF_THEN_ELSE(Time > 2020, 234.638, 234.638),
            0 ~ (3000000 + Time_less_Greenland_slide_experiment_start_yr) - Time,
            0 ~ Time_to_degrade_Kyoto_Flour_yr - ESCIMO_IF_THEN_ELSE(Time > 2020, 50.0, 50.0),
            0 ~ Time_to_regrow_NF_after_buning_yr - ESCIMO_IF_THEN_ELSE(Time > 2020, 30.0, 30.0),
            0 ~
              Tipping_point_search_emissions_GtCO2e_yr - ESCIMO_IF_THEN_ELSE(
                (Time >= 500000.0) & (Time <= Tipping_point_year_of_peak),
                Tipping_point_search_amount_at_start + (-Tipping_point_search_amount_at_start * (Time - 500000.0)) / (Tipping_point_year_of_peak - 500000.0),
                ESCIMO_IF_THEN_ELSE(Time > Tipping_point_year_of_peak, 0.0, InputEmissions_for_tipping_point_search),
              ),
            0 ~ Tipping_point_year_of_peak - 500001.0,
            0 ~
              (
                (((Total_carbon_in_Ocean_1850_GtC - Carbon_in_cold_ocean_0_to_100m_1850_GtC) - Carbon_in_cold_ocean_trunk_100m_to_bottom_1850_GtC) - Carbon_in_ocean_deep_1k_to_bottom_ocean_1850_GtC) -
                Carbon_in_ocean_upwelling_100m_to_1km_1850_GtC
              ) - Carbon_in_warm_ocean_0_to_100m_1850_GtC,
            0 ~
              (
                (((Total_carbon_in_ocean_GtC - C_in_cold_surface_water_GtC) - C_in_cold_water_trunk_downwelling_GtC) - C_in_deep_water_volume_1km_to_bottom_GtC) -
                C_in_intermediate_upwelling_water_100m_to_1km_GtC
              ) - C_in_warm_surface_water_GtC,
            0 ~ Total_CO2e_emissions_as_f_peak__GtCO2e_yr - Tipping_point_search_emissions_GtCO2e_yr,
            0 ~ Total_net_aerosol_forcing_ZJ_yr - Total_net_aerosol_forcings_W_m2 * UNIT_conversion_W_m2_earth_to_ZJ_yr,
            0 ~ (Total_net_aerosol_forcings_W_m2 - Anthropogenic_aerosol_forcing) - Model_Volcanic_aerosol_forcing_W_m2,
            0 ~ Total_sea_level_change_from_thermal_expansion_m - 1000.0 * (Sea_level_change_from_thermal_expansion_deep_m + Sea_level_change_from_thermal_expansion_surface_m),
            0 ~
              (
                (((Total_volume_of_ocean_water_GcubicM - Cold_surface_water_volume_Gm3) - Cold_water_volume_downwelling_Gm3) - Deep_water_volume_1km_to_4km_Gm3) -
                Intermediate_upwelling_water_volume_100m_to_1km_Gm3
              ) - Warm_surface_water_volume_Gm3,
            0 ~ TROP_being_deforested_Mkm2_yr - Stop_of_human_deforestation * TROP_deforestation_cutoff * TROP_historical_deforestation_pct_yr * TROP_with_normal_cover_Mkm2,
            0 ~ TROP_being_harvested_by_clear_cutting_Mkm2_yr - 0.5TROP_being_harvested_Mkm2_yr,
            0 ~ TROP_being_harvested_Mkm2_yr + (-1000.0TROP_Use_of_NF_biomass_for_energy_GtBiomass_yr) / TROP_living_biomass_densitiy_tBiomass_pr_km2,
            0 ~ TROP_being_harvested_normally_Mkm2_yr - 0.5TROP_being_harvested_Mkm2_yr,
            0 ~ TROP_Biomass_in_construction_material_being_burnt_GtBiomass_yr - 0.025TROP_Biomass_locked_in_construction_material_GtBiomass,
            0 ~ TROP_Biomass_in_construction_material_left_to_rot_GtBiomass_yr - 0.025TROP_Biomass_locked_in_construction_material_GtBiomass,
            0 ~ TROP_biomass_new_growing_GtBiomass___yr - 0.3333333333333333TROP_potential_less_actual_living_biomass_GtBiomass,
            0 ~ TROP_burning_Mkm2_yr - (0.003 * Effect_of_temperature_on_fire_incidence_dimensionless) * TROP_with_normal_cover_Mkm2,
            0 ~ TROP_Dead_biomass_decomposing_GtBiomass_yr + -TROP_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass / TROP_Time_to_decompose_undisturbed_dead_biomass_yr,
            0 ~ TROP_DeadB_and_SOM_tB_per_km2 - (8500.0 * Effect_of_CO2_on_new_biomass_growth) * Effect_of_temperature_on_new_biomass_growth_dimensionless,
            0 ~ TROP_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr - (0.001TROP_DeadB_and_SOM_tB_per_km2) * TROP_being_deforested_Mkm2_yr,
            0 ~ TROP_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr - (0.0001TROP_DeadB_and_SOM_tB_per_km2) * TROP_being_harvested_normally_Mkm2_yr,
            0 ~ TROP_deforestation_cutoff - var"combi_TROP_deforestation_cutoff_y[1]",
            0 ~ TROP_deforestation_cutoff_effect - var"combi_TROP_deforestation_cutoff_effect_y[1]",
            0 ~ TROP_deforested_as_pct_of_potial_area + -TROP_area_deforested_Mkm2 / TROP_potential_area_Mkm2,
            0 ~ TROP_deforestion_multiplier_wrt_2000 - var"combi_TROP_deforestion_multiplier_wrt_2000_y[1]",
            0 ~ TROP_for_construction_use_GtBiomass_yr - Use_of_TROP_biomass_for_construction_GtBiomass_yr,
            0 ~
              TROP_historical_deforestation_pct_yr -
              (0.01TROP_deforestion_multiplier_wrt_2000) * ESCIMO_IF_THEN_ELSE(false, EXP_12c_stopping_TROP_deforestation_from_2015, 1.0),
          ]
        end
        function generateEquations17()
          println("#Equation generated:" * "50" * "in: " * "generateEquations17")
          [
            0 ~ TROP_land_taken_out_of_use_GtBiomass - (0.001TROP_land_taken_out_of_use_Mkm2) * TROP_living_biomass_densitiy_tBiomass_pr_km2,
            0 ~ (((TROP_land_taken_out_of_use_Mkm2 - TROP_area_burnt_Mkm2) - TROP_area_clear_cut_Mkm2) - TROP_area_deforested_Mkm2) - TROP_area_harvested_Mkm2,
            0 ~ TROP_living_biomass_densitiy_tBiomass_pr_km2 - (16500.0 * Effect_of_CO2_on_new_biomass_growth) * Effect_of_temperature_on_new_biomass_growth_dimensionless,
            0 ~ TROP_Living_biomass_rotting_GtBiomass_yr - 0.016666666666666666TROP_Living_biomass_GtBiomass,
            0 ~
              TROP_NF_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr -
              (0.001TROP_living_biomass_densitiy_tBiomass_pr_km2) *
              (TROP_being_deforested_Mkm2_yr + TROP_being_harvested_by_clear_cutting_Mkm2_yr + TROP_being_harvested_normally_Mkm2_yr + TROP_burning_Mkm2_yr),
            0 ~ TROP_NF_regrowing_after_being_burnt_Mkm2_yr - 0.03333333333333333TROP_area_burnt_Mkm2,
            0 ~ TROP_NF_regrowing_after_harvesting_Mkm2_yr - 0.03333333333333333TROP_area_harvested_Mkm2,
            0 ~ (TROP_Living_biomass_GtBiomass + TROP_potential_less_actual_living_biomass_GtBiomass) - TROP_potential_living_biomass_GtBiomass,
            0 ~ TROP_potential_living_biomass_GtBiomass - (0.001TROP_living_biomass_densitiy_tBiomass_pr_km2) * (TROP_potential_area_Mkm2 - TROP_area_deforested_Mkm2),
            0 ~ TROP_regrowing_after_being_clear_cut_Mkm2_yr - 0.016666666666666666TROP_area_clear_cut_Mkm2,
            0 ~ TROP_regrowing_after_being_deforested_Mkm2_yr + -TROP_area_deforested_Mkm2 / Effective_Time_to_regrow_TROP_after_deforesting_yr,
            0 ~ TROP_runoff + -TROP_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass / TROP_runoff_time,
            0 ~ TROP_runoff_time - ESCIMO_IF_THEN_ELSE(Time > 2020, 2000.0, 2000.0),
            0 ~ TROP_soil_degradation_from_clear_cutting_GtBiomass_yr - (0.0005TROP_DeadB_and_SOM_tB_per_km2) * TROP_being_harvested_by_clear_cutting_Mkm2_yr,
            0 ~ TROP_soil_degradation_from_forest_fires_GtBiomass_yr,
            0 ~
              (
                (
                  (
                    (
                      (TROP_Sum_outflows_GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass_yr - TROP_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr) -
                      TROP_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr
                    ) - TROP_Dead_biomass_decomposing_GtBiomass_yr
                  ) - TROP_runoff
                ) - TROP_soil_degradation_from_clear_cutting_GtBiomass_yr
              ) - TROP_soil_degradation_from_forest_fires_GtBiomass_yr,
            0 ~ TROP_Time_to_decompose_undisturbed_dead_biomass_yr - ESCIMO_IF_THEN_ELSE(Time > 2020, 24.0, 24.0),
            0 ~ TROP_Use_of_NF_biomass_for_energy_GtBiomass_yr - Effect_of_population_and_urbanization_on_biomass_use * Use_of_TROP_for_energy_in_2000_GtBiomass,
            0 ~ (TROP_area_burnt_Mkm2 + TROP_area_clear_cut_Mkm2 + TROP_area_deforested_Mkm2 + TROP_area_harvested_Mkm2 + TROP_with_normal_cover_Mkm2) - TROP_potential_area_Mkm2,
            0 ~ TUNDRA_being_deforested_Mkm2_yr - Fraction_TUNDRA_being_deforested_1_yr * TUNDRA_with_normal_cover_Mkm2,
            0 ~ TUNDRA_being_harvested_Mkm2_yr + (-1000.0Use_of_TUNDRA_biomass_for_energy_GtBiomass_yr) / TUNDRA_living_biomass_densitiy_tBiomass_pr_km2,
            0 ~
              TUNDRA_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr -
              (0.001TUNDRA_living_biomass_densitiy_tBiomass_pr_km2) * (TUNDRA_being_deforested_Mkm2_yr + TUNDRA_being_harvested_Mkm2_yr + TUNDRA_burning_Mkm2_yr),
            0 ~ TUNDRA_Biomass_in_construction_material_being_burnt_GtBiomass_yr - 0.05TUNDRA_Biomass_locked_in_construction_material_GtBiomass,
            0 ~ TUNDRA_Biomass_in_construction_material_left_to_rot_GtBiomass_yr - 0.05TUNDRA_Biomass_locked_in_construction_material_GtBiomass,
            0 ~ TUNDRA_biomass_new_growing_GtBiomass___yr - 0.3333333333333333TUNDRA_potential_less_actual_living_biomass_GtBiomass,
            0 ~ TUNDRA_burning_Mkm2_yr - (0.01 * Effect_of_temperature_on_fire_incidence_dimensionless) * TUNDRA_with_normal_cover_Mkm2,
            0 ~ TUNDRA_Dead_biomass_decomposing_GtBiomass_yr - 0.001TUNDRA_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass,
            0 ~ TUNDRA_DeadB_and_SOM_tB_per_km2 - (65000.0 * Effect_of_CO2_on_new_biomass_growth) * Effect_of_temperature_on_new_biomass_growth_dimensionless,
            0 ~ TUNDRA_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr - (0.001TUNDRA_DeadB_and_SOM_tB_per_km2) * TUNDRA_being_deforested_Mkm2_yr,
            0 ~ TUNDRA_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr - (0.0001TUNDRA_DeadB_and_SOM_tB_per_km2) * TUNDRA_being_harvested_Mkm2_yr,
            0 ~ TUNDRA_for_construction_use_GtBiomass_yr - Use_of_TUNDRA_biomass_for_construction_GtBiomass_yr,
            0 ~ TUNDRA_historical_deforestation_pct_yr,
            0 ~ TUNDRA_land_taken_out_of_use_GtBiomass - (0.001TUNDRA_land_taken_out_of_use_Mkm2) * TUNDRA_living_biomass_densitiy_tBiomass_pr_km2,
            0 ~ ((TUNDRA_land_taken_out_of_use_Mkm2 - TUNDRA_area_burnt_Mkm2) - TUNDRA_area_harvested_Mkm2) - TUNDRA_deforested_Mkm2,
            0 ~ TUNDRA_living_biomass_densitiy_tBiomass_pr_km2 - (14500.0 * Effect_of_CO2_on_new_biomass_growth) * Effect_of_temperature_on_new_biomass_growth_dimensionless,
            0 ~ TUNDRA_Living_biomass_rotting_GtBiomass_yr - 0.01TUNDRA_Living_biomass_GtBiomass,
            0 ~ (TUNDRA_Living_biomass_GtBiomass + TUNDRA_potential_less_actual_living_biomass_GtBiomass) - TUNDRA_potential_living_biomass_GtBiomass,
            0 ~ TUNDRA_potential_living_biomass_GtBiomass - (0.001TUNDRA_living_biomass_densitiy_tBiomass_pr_km2) * (Tundra_potential_area_Mkm2 - TUNDRA_deforested_Mkm2),
            0 ~ TUNDRA_regrowing_after_being_burnt_Mkm2_yr - 0.1TUNDRA_area_burnt_Mkm2,
            0 ~ TUNDRA_regrowing_after_being_deforested_Mkm2_yr - 0.0125TUNDRA_deforested_Mkm2,
            0 ~ TUNDRA_regrowing_after_harvesting_Mkm2_yr - 0.1TUNDRA_area_harvested_Mkm2,
            0 ~ TUNDRA_runoff - 0.0005TUNDRA_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass,
            0 ~ TUNDRA_soil_degradation_from_forest_fires_GtBiomass_yr,
            0 ~
              (
                (
                  (
                    (TUNDRA_Sum_outflows_GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass_yr - TUNDRA_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr) -
                    TUNDRA_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr
                  ) - TUNDRA_Dead_biomass_decomposing_GtBiomass_yr
                ) - TUNDRA_runoff
              ) - TUNDRA_soil_degradation_from_forest_fires_GtBiomass_yr,
            0 ~ (TUNDRA_area_burnt_Mkm2 + TUNDRA_area_harvested_Mkm2 + TUNDRA_deforested_Mkm2 + TUNDRA_with_normal_cover_Mkm2) - Tundra_potential_area_Mkm2,
            0 ~ UNIT_conversion_for_CH4_from_CO2e_to_C - 0.030000000000000006,
            0 ~ UNIT_conversion_for_CO2_from_CO2e_to_C - 0.2727272727272727,
            0 ~ UNIT_conversion_from_MtCH4_to_GtC - 0.00075,
            0 ~ UNIT_conversion_GtCO2e_to_GtC - 0.2727272727272727,
            0 ~ UNIT_conversion_mm_to_m - 0.001,
          ]
        end
        function generateEquations18()
          println("#Equation generated:" * "50" * "in: " * "generateEquations18")
          [
            0 ~ UNIT_conversion_W_m2_earth_to_ZJ_yr - 5.1e-7Seconds_per_yr,
            0 ~ UNIT_converter_GtC_Gm3_to_ymoles_litre - 8.326394671107411e7,
            0 ~ (4.0 + Upper_to_deep_ocean_temp_diff_in_1850_degC) - Temp__ocean__surface_in_1850_C,
            0 ~ Upwelling_from_deep - Flow_of_cold_water_sinking_to_very_bottom_GcubicM_per_yr,
            0 ~ Upwelling_to_surface - Flow_of_cold_water_sinking_to_very_bottom_GcubicM_per_yr,
            0 ~ Urban_area_fraction - max(0.0, min(1.0, 0.0006557377049180329var"combi_Population_Lookup_bn_y[1]")),
            0 ~ Urban_Mkm2 - (0.30000000000000004Area_of_earth_Mkm2) * Urban_area_fraction,
            0 ~ Urbanzation_Effect_on_biomass_use - var"combi_Urbanzation_Effect_on_biomass_use_y[1]",
            0 ~ Use_of_GRASS_biomass_for_construction_GtBiomass_yr - Effect_of_population_and_urbanization_on_biomass_use * Use_of_GRASS_for_construction_in_2000_GtBiomass,
            0 ~ Use_of_GRASS_biomass_for_energy_GtBiomass_yr - Effect_of_population_and_urbanization_on_biomass_use * Use_of_GRASS_for_energy_in_2000_GtBiomass,
            0 ~ Use_of_GRASS_for_construction_in_2000_GtBiomass - 0.155,
            0 ~ Use_of_GRASS_for_energy_in_2000_GtBiomass - 3.1,
            0 ~ Use_of_NF_biomass_for_construction_GtBiomass_yr - Effect_of_population_and_urbanization_on_biomass_use * Use_of_NF_for_construction_in_2000_GtBiomass,
            0 ~ Use_of_NF_biomass_for_energy_GtBiomass_yr - Effect_of_population_and_urbanization_on_biomass_use * Use_of_NF_for_energy_in_2000_GtBiomass,
            0 ~ Use_of_NF_for_construction_in_2000_GtBiomass - 0.6669999999999999,
            0 ~ Use_of_NF_for_energy_in_2000_GtBiomass - 1.2535,
            0 ~ Use_of_TROP_biomass_for_construction_GtBiomass_yr - Effect_of_population_and_urbanization_on_biomass_use * Use_of_TROP_for_construction_in_2000_GtBiomass,
            0 ~ Use_of_TROP_for_construction_in_2000_GtBiomass - 1.776,
            0 ~ Use_of_TROP_for_energy_in_2000_GtBiomass - 0.259,
            0 ~ Use_of_TUNDRA_biomass_for_construction_GtBiomass_yr - Effect_of_population_and_urbanization_on_biomass_use * Use_of_TUNDRA_for_construction_in_2000_GtBiomass,
            0 ~ Use_of_TUNDRA_biomass_for_energy_GtBiomass_yr - Effect_of_population_and_urbanization_on_biomass_use * Use_of_TUNDRA_for_energy_in_2000_GtBiomass,
            0 ~ Use_of_TUNDRA_for_construction_in_2000_GtBiomass - 0.15,
            0 ~ Use_of_TUNDRA_for_energy_in_2000_GtBiomass - 3.0,
            0 ~
              Volcanic_aerosols_emissions -
              ESCIMO_IF_THEN_ELSE(false, 0.0, ESCIMO_IF_THEN_ELSE(Time < 2008, -Historical_aerosol_forcing_volcanic, 0.0)),
            0 ~ Volcanic_aerosols_removed_from_stratosphere - Volcanic_aerosols_in_stratosphere,
            0 ~ Volume_cold_ocean_0_to_100m - 3.619e7Fraction_of_ocean_classified_as_cold_surface,
            0 ~ Volume_cold_ocean_downwelling_100m_to_bottom - 1.30284e9Fraction_of_ocean_classified_as_cold_surface,
            0 ~ Volume_expansion_from_thermal_expansion_deep_Gm3_km3 - (0.0099Deep_ocean__cold__volume) * (Thermal_expansion_deep_pct - Thermal_expansion_deep_in_1850_pct),
            0 ~ Volume_expansion_from_thermal_expansion_surface_Gm3_km3 - (0.009980000000000001Surface_ocean__warm__volume) * (Thermal_expansion_surface_pct - Thermal_expansion_surface_in_1850_pct),
            0 ~ Volume_ocean_deep_1km_to_bottom - 8.10656e8,
            0 ~ Volume_ocean_upwelling_100m_to_1km - 2.31616e8,
            0 ~
              ((((Volume_of_total_ocean_Gm3 - Volume_cold_ocean_0_to_100m) - Volume_cold_ocean_downwelling_100m_to_bottom) - Volume_ocean_deep_1km_to_bottom) - Volume_ocean_upwelling_100m_to_1km) -
              Volume_warm_ocean_0_to_100m,
            0 ~ Volume_warm_ocean_0_to_100m - 2.8952e7,
            0 ~ Warming_due_to_CH4_blocking_W_m2 - Blocked_by_CH4 * LW_TOA_radiation_from_atm_to_space_difference_wrt_1850_W_m2,
            0 ~ Warming_due_to_CO2_blocking_W_m2 - Blocked_by_CO2 * LW_TOA_radiation_from_atm_to_space_difference_wrt_1850_W_m2,
            0 ~ Warming_due_to_othGHG_blocking_W_m2 - Blocked_by_otherGHG * LW_TOA_radiation_from_atm_to_space_difference_wrt_1850_W_m2,
            0 ~ Warming_due_to_water_vapor_blocking_W_m2 - Blocked_by_H20 * LW_TOA_radiation_from_atm_to_space_difference_wrt_1850_W_m2,
            0 ~ Years_of_exponential_rise_dless - Years_of_exponential_rise_yr,
            0 ~ Years_of_exponential_rise_yr - ESCIMO_IF_THEN_ELSE(Time > 20000000, Time - 20000000, 0.0),
            0 ~ Years_still_needed_to_reach_zero_emission_goal_yr - 34.0,
            0 ~ (All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC + yr_on_yr_change_in_C_in_land_use_GtC_yr) - All_C_taken_out_due_to_change_in_land_use_GtC,
            0 ~ yr_on_yr_change_in_C_in_ocean_GtC_yr - ESCIMO_IF_THEN_ELSE(Time > 1860, Total_carbon_in_ocean_GtC - C_in_ocean_1_yr_ago_GtC, 0.0),
            0 ~ flow_NF_TUNDRA_Biomass_in_construction_material_left_to_rot_GtBiomass_yr - NF_TUNDRA_Biomass_in_construction_material_left_to_rot_GtBiomass_yr,
            0 ~ flow_NF_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr - NF_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr,
            0 ~ flow_Evaporation_aka_latent_heat_flow - Evaporation_aka_latent_heat_flow,
            0 ~ flow_C_runoff_from_biomass_soil - C_runoff_from_biomass_soil,
            0 ~ flow_Kyoto_Flour_degradation - Kyoto_Flour_degradation,
            0 ~ flow_N2O_degradation_MtN2O_yr - N2O_degradation_MtN2O_yr,
            0 ~ flow_LW_TOA_radiation_from_atm_to_space - LW_TOA_radiation_from_atm_to_space,
            0 ~ flow_TROP_Living_biomass_rotting_GtBiomass_yr - TROP_Living_biomass_rotting_GtBiomass_yr,
          ]
        end
        function generateEquations19()
          println("#Equation generated:" * "50" * "in: " * "generateEquations19")
          [
            0 ~ flow_CO2_flux_TUNDRA_to_atm_Gtc_yr - CO2_flux_TUNDRA_to_atm_Gtc_yr,
            0 ~ flow_Sifting_DESERT_to_GRASS_Mkm2_yr - Sifting_DESERT_to_GRASS_Mkm2_yr,
            0 ~ flow_Upwelling_from_deep - Upwelling_from_deep,
            0 ~ flow_TUNDRA_regrowing_after_being_burnt_Mkm2_yr - TUNDRA_regrowing_after_being_burnt_Mkm2_yr,
            0 ~ flow_TUNDRA_runoff - TUNDRA_runoff,
            0 ~ flow_Greenland_ice_that_slid_into_the_ocean_melting__pos__or_freezing__neg__km3_yr - Greenland_ice_that_slid_into_the_ocean_melting__pos__or_freezing__neg__km3_yr,
            0 ~ flow_NF_being_harvested_by_clear_cutting_Mkm2_yr - NF_being_harvested_by_clear_cutting_Mkm2_yr,
            0 ~
              flow_Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_arctic_ice_ZJ_yr -
              Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_arctic_ice_ZJ_yr,
            0 ~ flow_TROP_Biomass_in_construction_material_being_burnt_GtBiomass_yr - TROP_Biomass_in_construction_material_being_burnt_GtBiomass_yr,
            0 ~ flow_Glacial_ice_melting__pos__or_freezing__neg__km3_yr - Glacial_ice_melting__pos__or_freezing__neg__km3_yr,
            0 ~ flow_NATURE_CCS_Fig3_GtC_yr - NATURE_CCS_Fig3_GtC_yr,
            0 ~ flow_NF_biomass_new_growing_GtBiomass___yr - NF_biomass_new_growing_GtBiomass___yr,
            0 ~ flow_LW_clear_sky_emissions_to_surface - LW_clear_sky_emissions_to_surface,
            0 ~ flow_CH4_in_the_atmosphere_converted_to_CO2 - CH4_in_the_atmosphere_converted_to_CO2,
            0 ~ flow_NF_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr - NF_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr,
            0 ~ flow_TROP_biomass_new_growing_GtBiomass___yr - TROP_biomass_new_growing_GtBiomass___yr,
            0 ~ flow_GRASS_Living_biomass_rotting_GtBiomass_yr - GRASS_Living_biomass_rotting_GtBiomass_yr,
            0 ~ flow_TUNDRA_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr - TUNDRA_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr,
            0 ~ flow_TUNDRA_biomass_new_growing_GtBiomass___yr - TUNDRA_biomass_new_growing_GtBiomass___yr,
            0 ~ flow_CO2_flux_from_atm_to_NF_for_new_growth_GtC_yr - CO2_flux_from_atm_to_NF_for_new_growth_GtC_yr,
            0 ~ flow_CH4_conversion_to_CO2_and_H2O - CH4_conversion_to_CO2_and_H2O,
            0 ~ flow_Flow_of_heat_to_deep_ocean_btw_72_and_08 - Flow_of_heat_to_deep_ocean_btw_72_and_08,
            0 ~ flow_GRASS_for_construction_use_GtBiomass_yr - GRASS_for_construction_use_GtBiomass_yr,
            0 ~ flow_Flow_of_cold_surface_water_welling_down_GcubicM_per_yr - Flow_of_cold_surface_water_welling_down_GcubicM_per_yr,
            0 ~ flow_TROP_for_construction_use_GtBiomass_yr - TROP_for_construction_use_GtBiomass_yr,
            0 ~ flow_Annual_glacial_ice_losing__pos__or_gaining__neg__GtIce_yr - Annual_glacial_ice_losing__pos__or_gaining__neg__GtIce_yr,
            0 ~ flow_NF_runoff - NF_runoff,
            0 ~ flow_NF_soil_degradation_from_forest_fires_GtBiomass_yr - NF_soil_degradation_from_forest_fires_GtBiomass_yr,
            0 ~ flow_GRASS_runoff - GRASS_runoff,
            0 ~ flow_Greenland_ice_sliding_into_the_ocean_km3_yr - Greenland_ice_sliding_into_the_ocean_km3_yr,
            0 ~ flow_TUNDRA_soil_degradation_from_forest_fires_GtBiomass_yr - TUNDRA_soil_degradation_from_forest_fires_GtBiomass_yr,
            0 ~ flow_Greenland_ice_melting__pos__or_freezing__neg__km3_yr - Greenland_ice_melting__pos__or_freezing__neg__km3_yr,
            0 ~ flow_SW_surface_absorption - SW_surface_absorption,
            0 ~ flow_All_N2O_emissions_MtN2O_yr - All_N2O_emissions_MtN2O_yr,
            0 ~ flow_NF_being_harvested_normally_Mkm2_yr - NF_being_harvested_normally_Mkm2_yr,
            0 ~ flow_Kyoto_Flour_emissions - Kyoto_Flour_emissions,
            0 ~ flow_CO2_flux_from_atm_to_TUNDRA_for_new_growth_GtC_yr - CO2_flux_from_atm_to_TUNDRA_for_new_growth_GtC_yr,
            0 ~ flow_Shifting_NF_to_TROP_Mkm2_yr - Shifting_NF_to_TROP_Mkm2_yr,
            0 ~
              flow_GRASS_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr -
              GRASS_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr,
            0 ~ flow_TUNDRA_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr - TUNDRA_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr,
            0 ~ flow_Carbon_flow_from_warm_to_cold_surface_GtC_per_yr - Carbon_flow_from_warm_to_cold_surface_GtC_per_yr,
            0 ~ flow_Shifting_GRASS_to_DESERT_Mkm2_yr - Shifting_GRASS_to_DESERT_Mkm2_yr,
            0 ~ flow_NF_being_deforested_Mkm2_yr - NF_being_deforested_Mkm2_yr,
            0 ~ flow_Heat_withdrawn_from_atm_by_melting__pos____added__neg__by_freezing_ice_ZJ_yr - Heat_withdrawn_from_atm_by_melting__pos____added__neg__by_freezing_ice_ZJ_yr,
            0 ~ flow_GRASS_biomass_new_growing_GtBiomass___yr - GRASS_biomass_new_growing_GtBiomass___yr,
            0 ~ flow_Man_made_fossil_C_emissions_GtC_yr - Man_made_fossil_C_emissions_GtC_yr,
            0 ~ flow_Greenland_ice_melting_as_water_km3_yr - Greenland_ice_melting_as_water_km3_yr,
            0 ~ flow_TROP_runoff - TROP_runoff,
            0 ~ flow_Antarctic_ice_losing__pos__or_gaining__neg__GtIce_yr - Antarctic_ice_losing__pos__or_gaining__neg__GtIce_yr,
            0 ~ flow_Carbon_flow_from_cold_surface_downwelling_Gtc_per_yr - Carbon_flow_from_cold_surface_downwelling_Gtc_per_yr,
          ]
        end
        function generateEquations20()
          println("#Equation generated:" * "50" * "in: " * "generateEquations20")
          [
            0 ~ flow_NF_regrowing_after_harvesting_Mkm2_yr - NF_regrowing_after_harvesting_Mkm2_yr,
            0 ~ flow_TROP_Dead_biomass_decomposing_GtBiomass_yr - TROP_Dead_biomass_decomposing_GtBiomass_yr,
            0 ~ flow_TUNDRA_being_deforested_Mkm2_yr - TUNDRA_being_deforested_Mkm2_yr,
            0 ~ flow_Shifting_TROP_to_GRASS_Mkm2_yr - Shifting_TROP_to_GRASS_Mkm2_yr,
            0 ~ flow_CH4_release___capture_from_permafrost_area_loss___gain_GtC_yr - CH4_release___capture_from_permafrost_area_loss___gain_GtC_yr,
            0 ~ flow_Volcanic_aerosols_emissions - Volcanic_aerosols_emissions,
            0 ~ flow_GRASS_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr - GRASS_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr,
            0 ~ flow_Natural_CH4_emissions - Natural_CH4_emissions,
            0 ~ flow_Flow_of_heat_to_atm_ZJ_yr - Flow_of_heat_to_atm_ZJ_yr,
            0 ~ flow_NF_Biomass_in_construction_material_being_burnt_GtBiomass_yr - NF_Biomass_in_construction_material_being_burnt_GtBiomass_yr,
            0 ~ flow_Flow_of_heat_to_deep_ocean - Flow_of_heat_to_deep_ocean,
            0 ~ flow_LW_surface_emission - LW_surface_emission,
            0 ~ flow_NF_regrowing_after_being_burnt_Mkm2_yr - NF_regrowing_after_being_burnt_Mkm2_yr,
            0 ~
              flow_TROP_NF_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr -
              TROP_NF_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr,
            0 ~ flow_TUNDRA_Biomass_in_construction_material_being_burnt_GtBiomass_yr - TUNDRA_Biomass_in_construction_material_being_burnt_GtBiomass_yr,
            0 ~ flow_Man_made_fossil_C_emissions_for_cumulation_GtC_yr - Man_made_fossil_C_emissions_for_cumulation_GtC_yr,
            0 ~ flow_C_absorption_by_ocean_from_atm_for_accumulation - C_absorption_by_ocean_from_atm_for_accumulation,
            0 ~ flow_Antarctic_ice_melting__pos__or_freezing__neg__km3_yr - Antarctic_ice_melting__pos__or_freezing__neg__km3_yr,
            0 ~ flow_Annual_flux_of_C_to_biomass_GtC_pr_yr - Annual_flux_of_C_to_biomass_GtC_pr_yr,
            0 ~ flow_GRASS_Biomass_in_construction_material_being_burnt_GtBiomass_yr - GRASS_Biomass_in_construction_material_being_burnt_GtBiomass_yr,
            0 ~ flow_NF_regrowing_after_being_deforested_Mkm2_yr - NF_regrowing_after_being_deforested_Mkm2_yr,
            0 ~ flow_Greenland_ice_losing__pos__or_gaining__neg__GtIce_yr - Greenland_ice_losing__pos__or_gaining__neg__GtIce_yr,
            0 ~ flow_CO2_flux_from_atm_to_TROP_for_new_growth_GtC_yr - CO2_flux_from_atm_to_TROP_for_new_growth_GtC_yr,
            0 ~ flow_NF_soil_degradation_from_clear_cutting_GtBiomass_yr - NF_soil_degradation_from_clear_cutting_GtBiomass_yr,
            0 ~ flow_Annual_release_of_C_from_permafrost_GtC_y - Annual_release_of_C_from_permafrost_GtC_y,
            0 ~ flow_Avg_volcanic_activity_GtC_yr - Avg_volcanic_activity_GtC_yr,
            0 ~ flow_TUNDRA_regrowing_after_harvesting_Mkm2_yr - TUNDRA_regrowing_after_harvesting_Mkm2_yr,
            0 ~ flow_Shifting_ice_on_land_to_tundra_Mkm2_yr - Shifting_ice_on_land_to_tundra_Mkm2_yr,
            0 ~ flow_C_diffusion_into_ocean_from_atm - C_diffusion_into_ocean_from_atm,
            0 ~ flow_Glacial_ice_melting_as_water_km3_yr - Glacial_ice_melting_as_water_km3_yr,
            0 ~ flow_NF_for_construction_use_GtBiomass_yr - NF_for_construction_use_GtBiomass_yr,
            0 ~ flow_Flow_of_heat_to_surface_ocean - Flow_of_heat_to_surface_ocean,
            0 ~ flow_TROP_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr - TROP_DeadB_SOM_being_lost_due_to_energy_harvesting_GtBiomass_yr,
            0 ~ flow_C_released_from_permafrost_as_either_CH4_or_CO2_in_GtC_yr - C_released_from_permafrost_as_either_CH4_or_CO2_in_GtC_yr,
            0 ~ flow_TROP_soil_degradation_from_forest_fires_GtBiomass_yr - TROP_soil_degradation_from_forest_fires_GtBiomass_yr,
            0 ~ flow_TROP_being_harvested_by_clear_cutting_Mkm2_yr - TROP_being_harvested_by_clear_cutting_Mkm2_yr,
            0 ~ flow_NF_regrowing_after_being_clear_cut_Mkm2_yr - NF_regrowing_after_being_clear_cut_Mkm2_yr,
            0 ~ flow_GRASS_being_harvested_Mkm2_yr - GRASS_being_harvested_Mkm2_yr,
            0 ~ flow_Convection_aka_sensible_heat_flow - Convection_aka_sensible_heat_flow,
            0 ~ flow_TUNDRA_for_construction_use_GtBiomass_yr - TUNDRA_for_construction_use_GtBiomass_yr,
            0 ~ flow_NF_burning_Mkm2_yr - NF_burning_Mkm2_yr,
            0 ~
              flow_NF_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr -
              NF_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr,
            0 ~ flow_TUNDRA_burning_Mkm2_yr - TUNDRA_burning_Mkm2_yr,
            0 ~ flow_CO2_flux_TROP_to_atm_GtC_yr - CO2_flux_TROP_to_atm_GtC_yr,
            0 ~ flow_Shifting_tundra_to_ice_on_land_Mkm2_yr - Shifting_tundra_to_ice_on_land_Mkm2_yr,
            0 ~ flow_Flow_of_water_from_warm_to_cold_surface_Gcubicm_per_yr - Flow_of_water_from_warm_to_cold_surface_Gcubicm_per_yr,
            0 ~ flow_Shifting_Tundra_to_NF_Mkm2_yr - Shifting_Tundra_to_NF_Mkm2_yr,
            0 ~ flow_Flow_of_heat_to_surface_ocean_btw_1972_and_2008 - Flow_of_heat_to_surface_ocean_btw_1972_and_2008,
            0 ~ flow_TUNDRA_Living_biomass_rotting_GtBiomass_yr - TUNDRA_Living_biomass_rotting_GtBiomass_yr,
            0 ~ flow_Methanehydrate_experimental_release_GtC__yr - Methanehydrate_experimental_release_GtC__yr,
          ]
        end
        function generateEquations21()
          println("#Equation generated:" * "50" * "in: " * "generateEquations21")
          [
            0 ~ flow_GRASS_regrowing_after_being_burnt_Mkm2_yr - GRASS_regrowing_after_being_burnt_Mkm2_yr,
            0 ~ flow_Montreal_gases_degradation - Montreal_gases_degradation,
            0 ~ flow_Carbon_flow_from_cold_to_deep_GtC_per_yr - Carbon_flow_from_cold_to_deep_GtC_per_yr,
            0 ~ flow_GRASS_soil_degradation_from_forest_fires_GtBiomass_yr - GRASS_soil_degradation_from_forest_fires_GtBiomass_yr,
            0 ~ flow_Shifting_TROP_to_NF_Mkm2_yr - Shifting_TROP_to_NF_Mkm2_yr,
            0 ~ flow_GRASS_being_deforested_Mkm2_yr - GRASS_being_deforested_Mkm2_yr,
            0 ~ flow_Shifting_GRASS_to_NF_Mkm2_yr - Shifting_GRASS_to_NF_Mkm2_yr,
            0 ~ flow_TROP_being_deforested_Mkm2_yr - TROP_being_deforested_Mkm2_yr,
            0 ~ flow_Arctic_ice_melting__pos__or_freezing__neg__km2_yr - Arctic_ice_melting__pos__or_freezing__neg__km2_yr,
            0 ~ flow_CO2_flux_from_atm_to_GRASS_for_new_growth_GtC_yr - CO2_flux_from_atm_to_GRASS_for_new_growth_GtC_yr,
            0 ~ flow_GRASS_regrowing_after_being_deforested_Mkm2_yr - GRASS_regrowing_after_being_deforested_Mkm2_yr,
            0 ~ flow_Net_C_to_atm_rate - Net_C_to_atm_rate,
            0 ~ flow_Methane_hydrates_released_and_converted_to_CO2_by_bacteria_GtC - Methane_hydrates_released_and_converted_to_CO2_by_bacteria_GtC,
            0 ~ flow_LW_surface_emissions_NOT_escaping_through_atm_window - LW_surface_emissions_NOT_escaping_through_atm_window,
            0 ~ flow_Antarctic_ice_melting_as_water_km3_yr - Antarctic_ice_melting_as_water_km3_yr,
            0 ~ flow_TUNDRA_Biomass_in_construction_material_left_to_rot_GtBiomass_yr - TUNDRA_Biomass_in_construction_material_left_to_rot_GtBiomass_yr,
            0 ~ flow_TROP_NF_regrowing_after_harvesting_Mkm2_yr - TROP_NF_regrowing_after_harvesting_Mkm2_yr,
            0 ~ flow_TUNDRA_being_harvested_Mkm2_yr - TUNDRA_being_harvested_Mkm2_yr,
            0 ~ flow_Net_heat_flow_ocean_from_surface_to_deep__ZJ_yr_ - Net_heat_flow_ocean_from_surface_to_deep__ZJ_yr_,
            0 ~ flow_TROP_regrowing_after_being_clear_cut_Mkm2_yr - TROP_regrowing_after_being_clear_cut_Mkm2_yr,
            0 ~ flow_GRASS_Biomass_in_construction_material_left_to_rot_GtBiomass_yr - GRASS_Biomass_in_construction_material_left_to_rot_GtBiomass_yr,
            0 ~ flow_TROP_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr - TROP_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr,
            0 ~ flow_Carbon_flow_from_deep - Carbon_flow_from_deep,
            0 ~ flow_Rate_of_destruction_of_wetlands - Rate_of_destruction_of_wetlands,
            0 ~ flow_Montreal_gases_emissions - Montreal_gases_emissions,
            0 ~ flow_LW_re_radiated_by_clouds - LW_re_radiated_by_clouds,
            0 ~
              flow_Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_antarctic_ice_ZJ_yr -
              Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_antarctic_ice_ZJ_yr,
            0 ~ flow_Depositing_of_C_to_sediment - Depositing_of_C_to_sediment,
            0 ~ flow_TUNDRA_Dead_biomass_decomposing_GtBiomass_yr - TUNDRA_Dead_biomass_decomposing_GtBiomass_yr,
            0 ~ flow_TUNDRA_regrowing_after_being_deforested_Mkm2_yr - TUNDRA_regrowing_after_being_deforested_Mkm2_yr,
            0 ~ flow_TROP_burning_Mkm2_yr - TROP_burning_Mkm2_yr,
            0 ~ flow_TROP_NF_regrowing_after_being_burnt_Mkm2_yr - TROP_NF_regrowing_after_being_burnt_Mkm2_yr,
            0 ~ flow_SW_Atmospheric_absorption - SW_Atmospheric_absorption,
            0 ~ flow_GRASS_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr - GRASS_DeadB_SOM_being_lost_due_to_deforestation_GtBiomass_yr,
            0 ~ flow_GRASS_regrowing_after_harvesting_Mkm2_yr - GRASS_regrowing_after_harvesting_Mkm2_yr,
            0 ~ flow_TROP_being_harvested_normally_Mkm2_yr - TROP_being_harvested_normally_Mkm2_yr,
            0 ~ flow_C_release_from_permafrost_melting_as_CO2_GtC_yr - C_release_from_permafrost_melting_as_CO2_GtC_yr,
            0 ~ flow_Human_activity_CH4_emissions - Human_activity_CH4_emissions,
            0 ~ flow_GRASS_Dead_biomass_decomposing_GtBiomass_yr - GRASS_Dead_biomass_decomposing_GtBiomass_yr,
            0 ~ flow_TROP_Biomass_in_construction_material_left_to_rot_GtBiomass_yr - TROP_Biomass_in_construction_material_left_to_rot_GtBiomass_yr,
            0 ~ flow_TROP_soil_degradation_from_clear_cutting_GtBiomass_yr - TROP_soil_degradation_from_clear_cutting_GtBiomass_yr,
            0 ~
              flow_TUNDRA_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr -
              TUNDRA_biomass_being_lost_from_deforestation__fires__energy_harvesting_and_clear_cutting_GtBiomass_yr,
            0 ~ flow_Shifting_NF_to_GRASS_Mkm2_yr - Shifting_NF_to_GRASS_Mkm2_yr,
            0 ~ flow_Heat_flow_from_the_earths_core - Heat_flow_from_the_earths_core,
            0 ~
              flow_Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_Greenland_ice_that_slid_into_the_ocean_ZJ_yr -
              Heat_withdrawn_from__ocean__surface_by_melting__pos____added__neg__by_freezing_Greenland_ice_that_slid_into_the_ocean_ZJ_yr,
            0 ~ flow_TROP_regrowing_after_being_deforested_Mkm2_yr - TROP_regrowing_after_being_deforested_Mkm2_yr,
            0 ~ flow_C_removal_rate_from_atm_for_nature_May_2020_GtC_y - C_removal_rate_from_atm_for_nature_May_2020_GtC_y,
            0 ~ flow_GRASS_burning_Mkm2_yr - GRASS_burning_Mkm2_yr,
            0 ~ flow_CO2_flux_GRASS_to_atm_Gtc_yr - CO2_flux_GRASS_to_atm_Gtc_yr,
            0 ~ flow_Upwelling_to_surface - Upwelling_to_surface,
          ]
        end
        function generateEquations22()
          println("#Equation generated:" * "50" * "in: " * "generateEquations22")
          [
            0 ~ flow_NF_Dead_biomass_decomposing_GtBiomass_yr - NF_Dead_biomass_decomposing_GtBiomass_yr,
            0 ~ flow_Carbon_captured_and_stored_GtC___yr - Carbon_captured_and_stored_GtC___yr,
            0 ~ flow_Volcanic_aerosols_removed_from_stratosphere - Volcanic_aerosols_removed_from_stratosphere,
            0 ~ flow_Carbon_flow_from_intermediate_to_surface_box_GtC_per_yr - Carbon_flow_from_intermediate_to_surface_box_GtC_per_yr,
            0 ~ flow_Greenland_ice_melting_that_slid_into_the_ocean_km3_yr - Greenland_ice_melting_that_slid_into_the_ocean_km3_yr,
            0 ~ flow_Shifting_NF_to_Tundra_Mkm2_yr - Shifting_NF_to_Tundra_Mkm2_yr,
            0 ~ flow_Shifting_GRASS_to_TROP_Mkm2_yr - Shifting_GRASS_to_TROP_Mkm2_yr,
            0 ~ flow_NF_Living_biomass_rotting_GtBiomass_yr - NF_Living_biomass_rotting_GtBiomass_yr,
            0 ~ flow_CO2_flux_NF_to_atm_Gtc_yr - CO2_flux_NF_to_atm_Gtc_yr,
            0 ~ flow_Flow_of_cold_water_sinking_to_very_bottom_GcubicM_per_yr - Flow_of_cold_water_sinking_to_very_bottom_GcubicM_per_yr,
            0 ~ flow_Biological_removal_of_C_from_WSW_GtC_per_yr - Biological_removal_of_C_from_WSW_GtC_per_yr,
            0 ~ Emissions_of_anthro_CH4_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_MtCH4_yr - var"combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_MtCH4_yr_y[1]",
            0 ~ combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_MtCH4_yr_u - Time,
            0 ~ Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP3_MtCH4_yr - var"combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP3_MtCH4_yr_y[1]",
            0 ~ combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP3_MtCH4_yr_u - Time,
            0 ~ Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP45_MtCH4_yr - var"combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP45_MtCH4_yr_y[1]",
            0 ~ combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP45_MtCH4_yr_u - Time,
            0 ~ Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP6_MtCH4_yr - var"combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP6_MtCH4_yr_y[1]",
            0 ~ combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP6_MtCH4_yr_u - Time,
            0 ~ Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP85_MtCH4_yr - var"combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP85_MtCH4_yr_y[1]",
            0 ~ combi_Emissions_of_anthro_CH4_from_Excel_1850_to_2100_RCP85_MtCH4_yr_u - Time,
            0 ~ Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_GtC_yr - var"combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_GtC_yr_y[1]",
            0 ~ combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_with_JR_2052_shape_GtC_yr_u - Time,
            0 ~
              Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_GtC_yr -
              var"combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_GtC_yr_y[1]",
            0 ~ combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_GtC_yr_u - Time,
            0 ~
              Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_hist_with_JR__large_scale_CCS__exp_GtC_yr -
              var"combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_hist_with_JR__large_scale_CCS__exp_GtC_yr_y[1]",
            0 ~ combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2015_RCP_hist_with_JR__large_scale_CCS__exp_GtC_yr_u - Time,
            0 ~ Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP3_GtC_yr - var"combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP3_GtC_yr_y[1]",
            0 ~ combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP3_GtC_yr_u - Time,
            0 ~ Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP45_GtC_yr - var"combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP45_GtC_yr_y[1]",
            0 ~ combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP45_GtC_yr_u - Time,
            0 ~ Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP6_GtC_yr - var"combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP6_GtC_yr_y[1]",
            0 ~ combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP6_GtC_yr_u - Time,
            0 ~ Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP85_GtC_yr - var"combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP85_GtC_yr_y[1]",
            0 ~ combi_Emissions_of_anthro_CO2_from_Excel_1850_to_2100_RCP85_GtC_yr_u - Time,
            0 ~
              Emissions_of_anthro_N20_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_MtN2O_yr -
              var"combi_Emissions_of_anthro_N20_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_MtN2O_yr_y[1]",
            0 ~ combi_Emissions_of_anthro_N20_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_MtN2O_yr_u - Time,
            0 ~ Emissions_of_anthro_N20_with_JR_2052_shape_MtN20_yr - var"combi_Emissions_of_anthro_N20_with_JR_2052_shape_MtN20_yr_y[1]",
            0 ~ combi_Emissions_of_anthro_N20_with_JR_2052_shape_MtN20_yr_u - Time,
            0 ~
              Emissions_of_Kyoto_Flour_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_kt_yr -
              var"combi_Emissions_of_Kyoto_Flour_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_kt_yr_y[1]",
            0 ~ combi_Emissions_of_Kyoto_Flour_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_kt_yr_u - Time,
            0 ~ Emissions_of_Kyoto_Flour_with_JR_2052_shape_kt_yr - var"combi_Emissions_of_Kyoto_Flour_with_JR_2052_shape_kt_yr_y[1]",
            0 ~ combi_Emissions_of_Kyoto_Flour_with_JR_2052_shape_kt_yr_u - Time,
            0 ~
              Emissions_of_Montreal_gases_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_kt_yr -
              var"combi_Emissions_of_Montreal_gases_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_kt_yr_y[1]",
            0 ~ combi_Emissions_of_Montreal_gases_from_Excel_1850_to_2015_RCP_hist_with_JR__twice_as_fast__exp_kt_yr_u - Time,
            0 ~ Emissions_of_Montreal_gases_with_JR_2052_shape_kt_yr - var"combi_Emissions_of_Montreal_gases_with_JR_2052_shape_kt_yr_y[1]",
            0 ~ combi_Emissions_of_Montreal_gases_with_JR_2052_shape_kt_yr_u - Time,
            0 ~ CH4_emissions_from_CO2e_C_Roads - var"combi_CH4_emissions_from_CO2e_C_Roads_y[1]",
            0 ~ combi_CH4_emissions_from_CO2e_C_Roads_u - Time,
            0 ~ CH4_emissions_from_CO2e_CAT - var"combi_CH4_emissions_from_CO2e_CAT_y[1]",
          ]
        end
        function generateEquations23()
          println("#Equation generated:" * "50" * "in: " * "generateEquations23")
          [
            0 ~ combi_CH4_emissions_from_CO2e_CAT_u - Time,
            0 ~ CH4_emissions_pct_contribution_to_Total_CO2e - var"combi_CH4_emissions_pct_contribution_to_Total_CO2e_y[1]",
            0 ~ combi_CH4_emissions_pct_contribution_to_Total_CO2e_u - Time,
            0 ~ CO2_emissions_from_CO2e_C_Roads - var"combi_CO2_emissions_from_CO2e_C_Roads_y[1]",
            0 ~ combi_CO2_emissions_from_CO2e_C_Roads_u - Time,
            0 ~ CO2_emissions_from_CO2e_CAT - var"combi_CO2_emissions_from_CO2e_CAT_y[1]",
            0 ~ combi_CO2_emissions_from_CO2e_CAT_u - Time,
            0 ~ CO2_emissions_pct_contribution_to_Total_CO2e - var"combi_CO2_emissions_pct_contribution_to_Total_CO2e_y[1]",
            0 ~ combi_CO2_emissions_pct_contribution_to_Total_CO2e_u - Time,
            0 ~ Historical_aerosol_emissions_anthro - var"combi_Historical_aerosol_emissions_anthro_y[1]",
            0 ~ combi_Historical_aerosol_emissions_anthro_u - Time,
            0 ~ Historical_forcing_from_solar_insolation_W_m2 - var"combi_Historical_forcing_from_solar_insolation_W_m2_y[1]",
            0 ~ combi_Historical_forcing_from_solar_insolation_W_m2_u - Time,
            0 ~ Historical_aerosol_forcing_volcanic - var"combi_Historical_aerosol_forcing_volcanic_y[1]",
            0 ~ combi_Historical_aerosol_forcing_volcanic_u - Time,
            0 ~ OGHG_Kyoto_Flour_emi_rcp3 - var"combi_OGHG_Kyoto_Flour_emi_rcp3_y[1]",
            0 ~ combi_OGHG_Kyoto_Flour_emi_rcp3_u - Time,
            0 ~ OGHG_Kyoto_Flour_emi_rcp45 - var"combi_OGHG_Kyoto_Flour_emi_rcp45_y[1]",
            0 ~ combi_OGHG_Kyoto_Flour_emi_rcp45_u - Time,
            0 ~ OGHG_Kyoto_Flour_emi_rcp6 - var"combi_OGHG_Kyoto_Flour_emi_rcp6_y[1]",
            0 ~ combi_OGHG_Kyoto_Flour_emi_rcp6_u - Time,
            0 ~ OGHG_Kyoto_Flour_emi_rcp85 - var"combi_OGHG_Kyoto_Flour_emi_rcp85_y[1]",
            0 ~ combi_OGHG_Kyoto_Flour_emi_rcp85_u - Time,
            0 ~ Kyoto_Flour_emissions_from_CO2e_C_Roads - var"combi_Kyoto_Flour_emissions_from_CO2e_C_Roads_y[1]",
            0 ~ combi_Kyoto_Flour_emissions_from_CO2e_C_Roads_u - Time,
            0 ~ Kyoto_Flour_emissions_from_CO2e_CAT - var"combi_Kyoto_Flour_emissions_from_CO2e_CAT_y[1]",
            0 ~ combi_Kyoto_Flour_emissions_from_CO2e_CAT_u - Time,
            0 ~ Kyoto_Flour_emissions_pct_contribution_to_Total_CO2e - var"combi_Kyoto_Flour_emissions_pct_contribution_to_Total_CO2e_y[1]",
            0 ~ combi_Kyoto_Flour_emissions_pct_contribution_to_Total_CO2e_u - Time,
            0 ~ OGHG_Montreal_gases_emi_rcp3 - var"combi_OGHG_Montreal_gases_emi_rcp3_y[1]",
            0 ~ combi_OGHG_Montreal_gases_emi_rcp3_u - Time,
            0 ~ OGHG_Montreal_gases_emi_rcp45 - var"combi_OGHG_Montreal_gases_emi_rcp45_y[1]",
            0 ~ combi_OGHG_Montreal_gases_emi_rcp45_u - Time,
            0 ~ OGHG_Montreal_gases_emi_rcp6 - var"combi_OGHG_Montreal_gases_emi_rcp6_y[1]",
            0 ~ combi_OGHG_Montreal_gases_emi_rcp6_u - Time,
            0 ~ OGHG_Montreal_gases_emi_rcp85 - var"combi_OGHG_Montreal_gases_emi_rcp85_y[1]",
            0 ~ combi_OGHG_Montreal_gases_emi_rcp85_u - Time,
            0 ~ othGHG_N20_man_made_emissions_rcp3 - var"combi_othGHG_N20_man_made_emissions_rcp3_y[1]",
            0 ~ combi_othGHG_N20_man_made_emissions_rcp3_u - Time,
            0 ~ othGHG_N20_man_made_emissions_rcp45 - var"combi_othGHG_N20_man_made_emissions_rcp45_y[1]",
            0 ~ combi_othGHG_N20_man_made_emissions_rcp45_u - Time,
            0 ~ othGHG_N20_man_made_emissions_rcp6 - var"combi_othGHG_N20_man_made_emissions_rcp6_y[1]",
            0 ~ combi_othGHG_N20_man_made_emissions_rcp6_u - Time,
            0 ~ othGHG_N20_man_made_emissions_rcp85 - var"combi_othGHG_N20_man_made_emissions_rcp85_y[1]",
            0 ~ combi_othGHG_N20_man_made_emissions_rcp85_u - Time,
            0 ~ RCP_3_CO2_concentration_1850_2100_ppm - var"combi_RCP_3_CO2_concentration_1850_2100_ppm_y[1]",
            0 ~ combi_RCP_3_CO2_concentration_1850_2100_ppm_u - Time,
            0 ~ RCP_45_CO2_concentration_1850_2100_ppm - var"combi_RCP_45_CO2_concentration_1850_2100_ppm_y[1]",
            0 ~ combi_RCP_45_CO2_concentration_1850_2100_ppm_u - Time,
            0 ~ RCP_6_CO2_concentration_1850_2100_ppm - var"combi_RCP_6_CO2_concentration_1850_2100_ppm_y[1]",
          ]
        end
        function generateEquations24()
          println("#Equation generated:" * "50" * "in: " * "generateEquations24")
          [
            0 ~ combi_RCP_6_CO2_concentration_1850_2100_ppm_u - Time,
            0 ~ RCP_85_CO2_concentration_1850_2100_ppm - var"combi_RCP_85_CO2_concentration_1850_2100_ppm_y[1]",
            0 ~ combi_RCP_85_CO2_concentration_1850_2100_ppm_u - Time,
            0 ~ Montreal_gases_emissions_from_CO2e_C_Roads - var"combi_Montreal_gases_emissions_from_CO2e_C_Roads_y[1]",
            0 ~ combi_Montreal_gases_emissions_from_CO2e_C_Roads_u - Time,
            0 ~ Montreal_gases_emissions_from_CO2e_CAT - var"combi_Montreal_gases_emissions_from_CO2e_CAT_y[1]",
            0 ~ combi_Montreal_gases_emissions_from_CO2e_CAT_u - Time,
            0 ~ Montreal_gases_emissions_pct_contribution_to_Total_CO2e - var"combi_Montreal_gases_emissions_pct_contribution_to_Total_CO2e_y[1]",
            0 ~ combi_Montreal_gases_emissions_pct_contribution_to_Total_CO2e_u - Time,
            0 ~ N2O_man_made_emissions_from_CO2e_C_Roads - var"combi_N2O_man_made_emissions_from_CO2e_C_Roads_y[1]",
            0 ~ combi_N2O_man_made_emissions_from_CO2e_C_Roads_u - Time,
            0 ~ N2O_man_made_emissions_from_CO2e_CAT - var"combi_N2O_man_made_emissions_from_CO2e_CAT_y[1]",
            0 ~ combi_N2O_man_made_emissions_from_CO2e_CAT_u - Time,
            0 ~ N2O_emissions_pct_contribution_to_Total_CO2e - var"combi_N2O_emissions_pct_contribution_to_Total_CO2e_y[1]",
            0 ~ combi_N2O_emissions_pct_contribution_to_Total_CO2e_u - Time,
            0 ~ Sea_level_rise_history_mm - var"combi_Sea_level_rise_history_mm_y[1]",
            0 ~ combi_Sea_level_rise_history_mm_u - Time,
            0 ~ E3_SC_1_CO2_GtC_yr - var"combi_E3_SC_1_CO2_GtC_yr_y[1]",
            0 ~ combi_E3_SC_1_CO2_GtC_yr_u - Time,
            0 ~ E3_SC_1_CH4_GtC_yr - var"combi_E3_SC_1_CH4_GtC_yr_y[1]",
            0 ~ combi_E3_SC_1_CH4_GtC_yr_u - Time,
            0 ~ E3_SC_1_N2O_Mt_yr - var"combi_E3_SC_1_N2O_Mt_yr_y[1]",
            0 ~ combi_E3_SC_1_N2O_Mt_yr_u - Time,
            0 ~ E3_SC_1_Kyoto_F_kt_yr - var"combi_E3_SC_1_Kyoto_F_kt_yr_y[1]",
            0 ~ combi_E3_SC_1_Kyoto_F_kt_yr_u - Time,
            0 ~ E3_SC_1_Montreal_gases_kt_yr - var"combi_E3_SC_1_Montreal_gases_kt_yr_y[1]",
            0 ~ combi_E3_SC_1_Montreal_gases_kt_yr_u - Time,
            0 ~ E3_SC_2_CO2_GtC_yr - var"combi_E3_SC_2_CO2_GtC_yr_y[1]",
            0 ~ combi_E3_SC_2_CO2_GtC_yr_u - Time,
            0 ~ E3_SC_2_CH4_GtC_yr - var"combi_E3_SC_2_CH4_GtC_yr_y[1]",
            0 ~ combi_E3_SC_2_CH4_GtC_yr_u - Time,
            0 ~ E3_SC_2_N2O_Mt_yr - var"combi_E3_SC_2_N2O_Mt_yr_y[1]",
            0 ~ combi_E3_SC_2_N2O_Mt_yr_u - Time,
            0 ~ E3_SC_2_Kyoto_F_kt_yr - var"combi_E3_SC_2_Kyoto_F_kt_yr_y[1]",
            0 ~ combi_E3_SC_2_Kyoto_F_kt_yr_u - Time,
            0 ~ E3_SC_2_Montreal_gases_kt_yr - var"combi_E3_SC_2_Montreal_gases_kt_yr_y[1]",
            0 ~ combi_E3_SC_2_Montreal_gases_kt_yr_u - Time,
            0 ~ E3_SC_3_CO2_GtC_yr - var"combi_E3_SC_3_CO2_GtC_yr_y[1]",
            0 ~ combi_E3_SC_3_CO2_GtC_yr_u - Time,
            0 ~ E3_SC_3_CH4_GtC_yr - var"combi_E3_SC_3_CH4_GtC_yr_y[1]",
            0 ~ combi_E3_SC_3_CH4_GtC_yr_u - Time,
            0 ~ E3_SC_3_N2O_Mt_yr - var"combi_E3_SC_3_N2O_Mt_yr_y[1]",
            0 ~ combi_E3_SC_3_N2O_Mt_yr_u - Time,
            0 ~ E3_SC_3_Kyoto_F_kt_yr - var"combi_E3_SC_3_Kyoto_F_kt_yr_y[1]",
            0 ~ combi_E3_SC_3_Kyoto_F_kt_yr_u - Time,
            0 ~ E3_SC_3_Montreal_gases_kt_yr - var"combi_E3_SC_3_Montreal_gases_kt_yr_y[1]",
            0 ~ combi_E3_SC_3_Montreal_gases_kt_yr_u - Time,
            0 ~ E3_SC_4_CO2_GtC_yr - var"combi_E3_SC_4_CO2_GtC_yr_y[1]",
            0 ~ combi_E3_SC_4_CO2_GtC_yr_u - Time,
            0 ~ E3_SC_4_CH4_GtC_yr - var"combi_E3_SC_4_CH4_GtC_yr_y[1]",
          ]
        end
        function generateEquations25()
          println("#Equation generated:" * "43" * "in: " * "generateEquations25")
          [
            0 ~ combi_E3_SC_4_CH4_GtC_yr_u - Time,
            0 ~ E3_SC_4_N2O_Mt_yr - var"combi_E3_SC_4_N2O_Mt_yr_y[1]",
            0 ~ combi_E3_SC_4_N2O_Mt_yr_u - Time,
            0 ~ E3_SC_4_Kyoto_F_kt_yr - var"combi_E3_SC_4_Kyoto_F_kt_yr_y[1]",
            0 ~ combi_E3_SC_4_Kyoto_F_kt_yr_u - Time,
            0 ~ E3_SC_4_Montreal_gases_kt_yr - var"combi_E3_SC_4_Montreal_gases_kt_yr_y[1]",
            0 ~ combi_E3_SC_4_Montreal_gases_kt_yr_u - Time,
            D(Model_N2O_concentration_in_1850_ppb) ~ 0.0,
            D(CO2_concentration_in_1850_ppm) ~ 0.0,
            D(Incoming_solar_in_1850_ZJ_yr) ~ 0.0,
            D(C_in_atmosphere_GtC_in_1850) ~ 0.0,
            D(C_in_biomass_in_1850_GtC) ~ 0.0,
            D(Total_carbon_in_ocean_GtC_in_1850) ~ 0.0,
            D(Temp_ocean_deep_1850_degC) ~ 0.0,
            D(init_ph_in_cold_water) ~ 0.0,
            D(Humidity_of_atmosphere_in_1850_g_kg) ~ 0.0,
            D(LW_TOA_radiation_from_atm_to_space_in_1850) ~ 0.0,
            D(Temp__ocean__surface_in_1850_C) ~ 0.0,
            D(Fraction_blocked_by_ALL_GHG_in_1850) ~ 0.0,
            D(Fraction_blocked_CO2_in_1850) ~ 0.0,
            D(Fraction_blocked_CH4_in_1850) ~ 0.0,
            D(Fraction_blocked_othGHG_in_1850) ~ 0.0,
            D(init_C_in_GRASS) ~ 0.0,
            D(init_C_in_NF) ~ 0.0,
            D(init_C_in_TROP) ~ 0.0,
            D(init_C_in_TUNDRA) ~ 0.0,
            D(Fossil_fuel_reserves_in_ground_1850_GtC) ~ 0.0,
#            D(Time) ~ 0.0,
            D(Aerosol_anthropogenic_emissions_in_2010) ~ 0.0,
            D(CO2_emissions_in_2010) ~ 0.0,
            D(CO2_ppm_value_at_When_to_sample) ~ 0.0,
            D(CO4_emissions_in_2010) ~ 0.0,
            D(Greenland_slide_experiment_end_condition) ~ 0.0,
            D(Kyoto_Flour_concentration_in_1970_ppt) ~ 0.0,
            D(Kyoto_Flour_emissions_RCPs_JR_in_2010) ~ 0.0,
            D(Montreal_gases_concentration_in_1970_ppt) ~ 0.0,
            D(Montreal_gases_emissions_RCPs_JR_in_2010) ~ 0.0,
            D(N20_emissions_RCPs_JR_in_2010) ~ 0.0,
            D(Tipping_point_search_amount_at_start) ~ 0.0,
            D(ifCond1) ~ 0.0,
            D(ifCond2) ~ 0.0,
            0 ~ ifelse(ifCond1 == true, 1.0 + 0.0003 * (Kyoto_Flour_concentration_ppt / Kyoto_Flour_concentration_in_1970_ppt - 1.0), 1.0) - ifEq_tmp304,
            0 ~ ifelse(ifCond2 == true, 1.0 + 0.003 * (Montreal_gases_concentration_ppt / Montreal_gases_concentration_in_1970_ppt - 1.0), 1.0) - ifEq_tmp305,
          ]
        end
        push!(equationConstructors, generateEquations0)
        push!(equationConstructors, generateEquations1)
        push!(equationConstructors, generateEquations2)
        push!(equationConstructors, generateEquations3)
        push!(equationConstructors, generateEquations4)
        push!(equationConstructors, generateEquations5)
        push!(equationConstructors, generateEquations6)
        push!(equationConstructors, generateEquations7)
        push!(equationConstructors, generateEquations8)
        push!(equationConstructors, generateEquations9)
        push!(equationConstructors, generateEquations10)
        push!(equationConstructors, generateEquations11)
        push!(equationConstructors, generateEquations12)
        push!(equationConstructors, generateEquations13)
        push!(equationConstructors, generateEquations14)
        push!(equationConstructors, generateEquations15)
        push!(equationConstructors, generateEquations16)
        push!(equationConstructors, generateEquations17)
        push!(equationConstructors, generateEquations18)
        push!(equationConstructors, generateEquations19)
        push!(equationConstructors, generateEquations20)
        push!(equationConstructors, generateEquations21)
        push!(equationConstructors, generateEquations22)
        push!(equationConstructors, generateEquations23)
        push!(equationConstructors, generateEquations24)
        push!(equationConstructors, generateEquations25)
      end
      for constructor in equationConstructors
        push!(equationComponents, constructor())
      end
      eqs = collect(Iterators.flatten(equationComponents))
      events = [Time > 1970 => [ifCond1 ~ true], !(Time > 1970) => [ifCond1 ~ false], Time > 1970 => [ifCond2 ~ true], !(Time > 1970) => [ifCond2 ~ false]]
      nonLinearSystem = ODESystem(eqs, t, vars, parameters; name = :($(Symbol("ESCIMO"))), discrete_events = events)
      firstOrderSystem = nonLinearSystem
      reducedSystem = firstOrderSystem
      local event_p = [
        0.0,
        0.7,
        0.7,
        0.17,
        0.7,
        0.24,
        0.4,
        0.4,
        0.08,
        0.3,
        0.16,
        0.7,
        0.13,
        0.18,
        0.08,
        0.1,
        0.168,
        0.14,
        0.23,
        0.23,
        0.23,
        0.15,
        0.0,
        4.0,
        0.0,
        0.0,
        3.0e7,
        0.7,
        1.34e7,
        15.0,
        0.2,
        0.4,
        17500.0,
        5.1e14,
        361900.0,
        0.0,
        0.0025,
        4.8e-5,
        0.1,
        1.0,
        2.14,
        2.14,
        1.35,
        600.0,
        1.69,
        0.5,
        2240.0,
        2240.0,
        2240.0,
        2240.0,
        2240.0,
        1720.81,
        7.3,
        35.0,
        0.2,
        0.071,
        468.0,
        0.04,
        0.04,
        1.0e-6,
        -1.325,
        2.8,
        -1.0,
        0.127044,
        0.916,
        5.0,
        0.19,
        1.0,
        1.0,
        0.289,
        EXP_12f_Stratospheric_scattering_experiment_0_off_1_on,
        5.0,
        0.0,
        30000.0,
        0.051,
        0.0837,
        0.006,
        0.158,
        1.0,
        1.0,
        0.7,
        0.6,
        0.5,
        0.1,
        0.9,
        0.8,
        167000.0,
        25.0,
        298.0,
        1.0,
        0.5,
        2.5,
        100.0,
        10.0,
        1.5,
        1200.0,
        0.5,
        1.0,
        0.1,
        0.0,
        14500.0,
        310.0,
        1.0,
        0.1,
        2000.0,
        2.0,
        1000.0,
        0.33,
        2.93e6,
        0.25,
        70.0,
        0.9167,
        0.0001717,
        1.9532e6,
        1025.67,
        25000.0,
        0.0003327,
        0.23,
        10.0,
        60.0,
        3.0,
        0.4,
        1.0,
        234.638,
        50.0,
        30.0,
        2000.0,
        24.0,
        273.15,
        7000.0,
        25.0,
        27.9,
        20.0,
        0.0398,
        0.303,
        10.0,
        35.0,
        0.71,
        10000.0,
        0.0594,
        5.35,
        0.12,
        363.504,
        900.0,
        9.0,
        0.4,
        NEvt_13a_double_rate_of_melting_ice_and_permafrost,
        NEvt_13b2_Double_incidence_of_biomass_fires,
        NEvt_13b3_double_sunspot_amplitude_from_2015_onwards_1_normal_2_double,
        NEvt_13c1_increase_in_area_covered_by_low_clouds,
        NEvt_13d_Greenland_slide_experiment_start_yr,
        NEvt_2a_Volcanic_eruptions_in_the_future_VAEs_first_future_pulse,
        NEvt_3b_increase_in_area_covered_by_high_clouds,
        2.5,
        0.0,
        1.0,
        20.0,
        3.0,
        330.0,
        27500.0,
        0.5,
        0.5,
        1.0,
        0.1,
        0.0,
        7500.0,
        115.0,
        0.7,
        0.02,
        2000.0,
        250.0,
        0.0,
        1.0,
        0.065,
        5.0,
        1.0,
        Policy_1_Reducing_GHG_emissions_by_one_third_by_2035,
        Policy_2_Large_scale_implementation_of_carbon_capture_and_geological_storage__CCS_,
        6.1,
        1.0,
        0.2,
        0.0,
        4.0,
        50.0,
        4.0,
        3.0,
        0.4,
        3.0,
        1.0,
        1.0,
        10.0,
        10000.0,
        273.15,
        0.23,
        0.220588,
        10.0,
        60.0,
        3.0,
        0.4,
        1.0,
        234.638,
        50.0,
        30.0,
        2000.0,
        24.0,
        1.0,
        2.5,
        0.58,
        50.0,
        50.0,
        58.0,
        5.0,
        0.0,
        0.0,
        0.0,
        0.3,
        0.3,
        0.1,
        1.0,
        1.0,
        2.0,
        0.1,
        1.0,
        5.0,
        0.1,
        0.2,
        0.01,
        0.2,
        0.05,
        0.2,
        5.0,
        0.1,
        1.2,
        0.65,
        0.1,
        0.71,
        0.1,
        0.05,
        -3.5,
        11.0,
        5.67037e-8,
        3.0e7,
        3.0,
        Switch_0_normal_model_1_dbl_CO2_2_1pct_incr,
        Switch_btw_historical_CO2_CH4_emissions_or_constant_1history_0constant,
        SWITCH_for_NATURE_comm_200115_base_1_cut_all_mm_emi_in_2020_2,
        SWITCH_future_slope_base_0_plus_5_1_minus_5_2,
        SWITCH_h2o_blocked_table_0_linear_1_poly_2,
        SWITCH_h2o_poly_dyn_0_equ_1,
        SWITCH_nature_rev_0_base_1_steeper_2_less_steep,
        Switch_to_choose_CO2_CH4_aerosol_other_GHG_0normal_1zero_from_2010_2constant_from_2010,
        Switch_to_choose_input_emission_scenario_for_CO2_CH4_and_oth_GHG,
        0.0,
        Switch_to_run_experiment_12a_reduction_in_emissions_0_off_1_on,
        Switch_to_run_experiment_12b_CCS_0_off_1_on,
        Switch_to_run_experiment_12c_stopping_TROP_deforestation_0_off_1_on,
        Switch_to_run_experiment_12e_white_surfaces_0_off_1_on,
        Switch_to_run_NATURE_experiment_CCS_0_off_1_on_0,
        Switch_to_run_POLICY_4_Stopping_logging_in_Northern_forests_0_off_1_on,
        4.0,
        274.31,
        9.7,
        286.815,
        2050.0,
        2800.0,
        800.0,
        100.0,
        3000.0,
        1.0,
        6.51772,
        739.89,
        211.397,
        26.227,
        30.0,
        95.0,
        20000.0,
        25.0,
        500.0,
        4000.0,
        20.0,
        18000.0,
        500.0,
        5.0,
        18.0,
        10.0,
        80.0,
        80.0,
        30.0,
        10.0,
        80.0,
        3.0,
        0.0,
        210000.0,
        500000.0,
        1.7,
        1.0,
        0.3,
        60.0,
        20.0,
        30.0,
        0.5,
        160.0,
        8500.0,
        0.5,
        0.5,
        1.0,
        0.1,
        0.0,
        16500.0,
        370.0,
        0.3,
        1.0,
        -0.5,
        3.0,
        2.0,
        0.0,
        2.5,
        100.0,
        10.0,
        1.5,
        1200.0,
        65000.0,
        0.5,
        1.0,
        0.1,
        0.0,
        14500.0,
        300.0,
        1.0,
        0.0,
        2000.0,
        3.0,
        1000.0,
        1.0,
        1.0,
        1.0,
        1.0e-6,
        1000.0,
        1.0,
        0.305,
        1.0,
        1.0e6,
        1000.0,
        1000.0,
        1000.0,
        1.0,
        1.0,
        1.0e6,
        1.0,
        1.0,
        1.0e6,
        1.0e12,
        31536.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        0.004,
        0.05,
        1.0,
        0.58,
        1.09,
        0.48,
        0.07,
        0.05,
        1.0,
        40.0,
        10.0,
        1.0,
        0.225,
        0.00125,
        1.0e7,
        When_first_destroyed_yr,
        When_methane_hydrates_first_released_yr,
        When_to_sample_for_CO2_experiment_yr,
        2020.0,
        273.15,
        1.0e21,
        2.0,
        2.0,
        2.0,
        2.0,
        2.0,
        2.0,
        2.0,
        2.0,
        2.0,
        2.0,
        2.0,
        2.0,
      ]
      local discreteVars = collect(
        values(
          ModelingToolkit.OrderedDict(
            Antarctic_ice_volume_km3 => 3.0e7,
            Arctic_ice__on_sea__area_km2 => 1.34e7,
            C_in_atmosphere_GtC => 600.0,
            C_in_atmosphere_in_form_of_CH4 => 1.69,
            C_in_cold_surface_water_GtC => Carbon_in_cold_ocean_0_to_100m_1850_GtC,
            C_in_cold_water_trunk_downwelling_GtC => Carbon_in_cold_ocean_trunk_100m_to_bottom_1850_GtC,
            C_in_deep_water_volume_1km_to_bottom_GtC => Carbon_in_ocean_deep_1k_to_bottom_ocean_1850_GtC,
            C_in_intermediate_upwelling_water_100m_to_1km_GtC => Carbon_in_ocean_upwelling_100m_to_1km_1850_GtC,
            C_in_permafrost_in_form_of_CH4 => 1200.0,
            C_in_sediment => 3.0e9,
            C_in_warm_surface_water_GtC => Carbon_in_warm_ocean_0_to_100m_1850_GtC,
            Cold_surface_water_volume_Gm3 => Volume_cold_ocean_0_to_100m,
            Cold_water_volume_downwelling_Gm3 => Volume_cold_ocean_downwelling_100m_to_bottom,
            Cumulative_antarctic_ice_volume_loss_GtIce => 0.0,
            Cumulative_C_released_from_permafrost_as_either_CH4_or_CO2_in_GtC_yr => 0.0,
            Cumulative_carbon_captured_and_stored_GtC => 0.0,
            Cumulative_carbon_removed_from_atm_for_nature_May_2020 => 0.0,
            Cumulative_flow_of_C_to_biomass_since_1850_GtC => 0.0,
            Cumulative_glacial_ice_volume_loss_GtIce => 0.0,
            Cumulative_Greenland_ice_volume_loss_GtIce => 0.0,
            Cumulative_heat_to_atm_ZJ => 0.0,
            Cumulative_ocean_volume_increase_due_to_ice_melting_km3 => 0.0,
            Cumulative_release_of_C_from_permafrost_GtC => 0.0,
            Deep_water_volume_1km_to_4km_Gm3 => Volume_ocean_deep_1km_to_bottom,
            DESERT_Mkm2 => 25.4,
            Fossil_fuel_reserves_in_ground_GtC => 6000.0,
            Glacial_ice_volume_km3 => 167000.0,
            GRASS_area_burnt_Mkm2 => 1.0,
            GRASS_area_harvested_Mkm2 => 2.5,
            GRASS_Biomass_locked_in_construction_material_GtBiomass => 1.5,
            GRASS_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass => 1200.0,
            GRASS_deforested_Mkm2 => 0.5,
            GRASS_Living_biomass_GtBiomass => 310.0,
            GRASS_potential_area_Mkm2 => 22.5,
            Greenland_ice_volume_on_Greenland_km3 => 2.93e6,
            Greenland_ice_volume_that_slid_into_the_ocean_km3 => 0.0,
            Heat_in_atmosphere_ZJ => 1025.67,
            Heat_in_deep_ZJ => 1.9532e6,
            Heat_in_surface => 25000.0,
            Intermediate_upwelling_water_volume_100m_to_1km_Gm3 => Volume_ocean_upwelling_100m_to_1km,
            Kyoto_Flour_gases_in_atm => 0.0,
            Montreal_gases_in_atm => 0.0,
            N2O_in_atmosphere_MtN2O => 900.0,
            NATURE_Cumulative_CCS_GtC => 0.0,
            NF_area_burnt_Mkm2 => 2.5,
            NF_area_clear_cut_Mkm2 => 1.0,
            NF_area_deforested_Mkm2 => 0.0,
            NF_area_harvested_Mkm2 => 1.0,
            NF_Biomass_locked_in_construction_material_GtBiomass => 3.0,
            NF_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass => 330.0,
            NF_Living_biomass_GtBiomass => 115.0,
            NF_potential_area_Mkm2 => 17.0,
            Sum_C_absorbed_by_ocean_GtC => 0.0,
            Sum_heat_to_deep_ocean => 0.0,
            Sum_heat_to_deep_ocean_btw_72_and_08 => 0.0,
            Sum_heat_to_surface_ocean_btw_72_and_08 => 0.0,
            Sum_heat_to_surface_ocean_ZJ => 0.0,
            Sum_man_made_CO2_emissions_GtC => 0.0,
            Sum_net_C_to_atm => 0.0,
            TROP_area_burnt_Mkm2 => 1.7,
            TROP_area_clear_cut_Mkm2 => 0.3,
            TROP_area_deforested_Mkm2 => 1.0,
            TROP_area_harvested_Mkm2 => 0.3,
            TROP_Biomass_locked_in_construction_material_GtBiomass => 30.0,
            TROP_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass => 160.0,
            TROP_Living_biomass_GtBiomass => 370.0,
            TROP_potential_area_Mkm2 => 25.0,
            TUNDRA_area_burnt_Mkm2 => 2.0,
            TUNDRA_area_harvested_Mkm2 => 2.5,
            TUNDRA_Biomass_locked_in_construction_material_GtBiomass => 1.5,
            TUNDRA_Dead_biomass__litter_and_soil_organic_matter_SOM_GtBiomass => 1200.0,
            TUNDRA_deforested_Mkm2 => 0.0,
            TUNDRA_Living_biomass_GtBiomass => 300.0,
            Tundra_potential_area_Mkm2 => 22.5,
            Volcanic_aerosols_in_stratosphere => 0.0,
            Warm_surface_water_volume_Gm3 => Volume_warm_ocean_0_to_100m,
            Wetlands_area => 1.0e7,
            Aerosol_anthropogenic_emissions_in_2010 => 0.0,
            CO2_emissions_in_2010 => 0.0,
            CO2_ppm_value_at_When_to_sample => MODEL_CO2_concentration_in_atmosphere2_ppm,
            CO4_emissions_in_2010 => 0.0,
            Greenland_slide_experiment_end_condition => 0.0,
            Kyoto_Flour_concentration_in_1970_ppt => 0.0,
            Kyoto_Flour_emissions_RCPs_JR_in_2010 => 0.0,
            Montreal_gases_concentration_in_1970_ppt => 0.0,
            Montreal_gases_emissions_RCPs_JR_in_2010 => 0.0,
            N20_emissions_RCPs_JR_in_2010 => 0.0,
            Tipping_point_search_amount_at_start => 12.0,
            Arctic_land_surface_temp_anomaly_compared_to_1850 => Temp_surface_anomaly_compared_to_1850_degC,
            Biological_removal_of_C_from_WSW_GtC_per_yr => Net_marine_primary_production_NMPP_GtC_pr_yr,
            Effect_of_temp_on_permafrost_melting_dmnl => 1.0 + Slope_btw_temp_and_permafrost_melting___freezing * (Temp_diff_relevant_for_melting_or_freezing_from_1850 / 4.0 - 1.0),
            Temp_diff_relevant_for_melting_or_freezing_arctic_ice_from_1850 => Temp_surface_anomaly_compared_to_1850_degC,
            Temp_diff_relevant_for_melting_or_freezing_from_1850 => Temp_surface_C - 13.66500000000002,
            yr_on_yr_change_in_C_in_atm_GtC_yr => C_in_atmosphere_GtC - C_in_atm_1_yr_ago_GtC,
            C_in_ocean_1_yr_ago_GtC => Total_carbon_in_ocean_GtC,
            C_in_ocean_1_yr_ago_GtC_LV1 => Total_carbon_in_ocean_GtC,
            C_in_ocean_1_yr_ago_GtC_LV2 => Total_carbon_in_ocean_GtC,
            Atmos_heat_used_for_melting_last_year_1_yr_LV => 0.0,
            Ocean_heat_used_for_melting_last_year_ZJ_yr_LV => 0.0,
            C_in_atm_1_yr_ago_GtC_LV3 => C_in_atm_1_yr_ago_GtC_DL * C_in_atmosphere_GtC,
            C_in_atm_1_yr_ago_GtC_LV2 => C_in_atm_1_yr_ago_GtC_LV3,
            C_in_atm_1_yr_ago_GtC_LV1 => C_in_atm_1_yr_ago_GtC_LV3,
            All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_LV3 => All_C_taken_out_due_to_change_in_land_use_GtC * All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_DL,
            All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_LV2 => All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_LV3,
            All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_LV1 => All_C_taken_out_due_to_change_in_land_use_GtC_1_yr_ago_GtC_LV3,
            Model_N2O_concentration_in_1850_ppb => 0.0,
            CO2_concentration_in_1850_ppm => 0.0,
            Incoming_solar_in_1850_ZJ_yr => 0.0,
            C_in_atmosphere_GtC_in_1850 => 0.0,
            C_in_biomass_in_1850_GtC => 0.0,
            Total_carbon_in_ocean_GtC_in_1850 => 0.0,
            Temp_ocean_deep_1850_degC => 0.0,
            init_ph_in_cold_water => 0.0,
            Humidity_of_atmosphere_in_1850_g_kg => 0.0,
            LW_TOA_radiation_from_atm_to_space_in_1850 => 0.0,
            Temp__ocean__surface_in_1850_C => 0.0,
            Fraction_blocked_by_ALL_GHG_in_1850 => 0.0,
            Fraction_blocked_CO2_in_1850 => 0.0,
            Fraction_blocked_CH4_in_1850 => 0.0,
            Fraction_blocked_othGHG_in_1850 => 0.0,
            init_C_in_GRASS => 0.0,
            init_C_in_NF => 0.0,
            init_C_in_TROP => 0.0,
            init_C_in_TUNDRA => 0.0,
            Fossil_fuel_reserves_in_ground_1850_GtC => 0.0,
            Time => 0.0,
            Aerosol_anthropogenic_emissions_in_2010 => 0.0,
            CO2_emissions_in_2010 => 0.0,
            CO2_ppm_value_at_When_to_sample => 0.0,
            CO4_emissions_in_2010 => 0.0,
            Greenland_slide_experiment_end_condition => 0.0,
            Kyoto_Flour_concentration_in_1970_ppt => 0.0,
            Kyoto_Flour_emissions_RCPs_JR_in_2010 => 0.0,
            Montreal_gases_concentration_in_1970_ppt => 0.0,
            Montreal_gases_emissions_RCPs_JR_in_2010 => 0.0,
            N20_emissions_RCPs_JR_in_2010 => 0.0,
            Tipping_point_search_amount_at_start => 0.0,
          ),
        ),
      )
      event_p = vcat(event_p, discreteVars)
      local aux = Vector{Any}(undef, 3)
      aux[1] = event_p
      aux[2] = Float64[]
      aux[3] = reducedSystem
      callbacks = ESCIMOCallbackSet(aux)
      problem = ModelingToolkit.ODEProblem(reducedSystem, initialValues, tspan, pars, callback = callbacks)
      return (problem, callbacks, initialValues, reducedSystem, tspan, pars, vars, irreductableSyms)
    end
  end
  function ESCIMOSimulate(tspan = (1850.0, 1900.0); solver = Rosenbrock23())
    (ESCIMOModel_problem, callbacks, ivs, ESCIMOModel_ReducedSystem, tspan, pars, vars, irreductable) = ESCIMOModel(tspan)
    solve(ESCIMOModel_problem, solver)
  end
