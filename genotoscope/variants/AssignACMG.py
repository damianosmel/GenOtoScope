import logging
from genotoscope.variants.VariantInfo import VariantInfo
from genotoscope.utils import split_name_assignment


class AssignACMG:
	"""
	Class to assign ACMG class for an input variant focusing on hearing loss.

	1) Possible class: pathogenic, likely pathogenic, benign, likely benign and uncertain significance
	Based on research work:
	Oza, Andrea M., et al. "Expert specification of the ACMG/AMP variant interpretation guidelines for
	genetic hearing loss." Human mutation 39.11 (2018): 1593-1613.
	DOI: 10.1002/humu.23630

	2) Compute variant pathogenicity probability based on the research work:
	Tavtigian, Sean V., et al. "Modeling the ACMG/AMP variant classification guidelines as a Bayesian
	classification framework." Genetics in Medicine 20.9 (2018): 1054-1060.
	DOI: 10.1038/gim.2017.210
	"""

	def __init__(self, variant_info, pathogenicity_odds_very_strong, scaling_factor, pathogenicity_prior):
		"""
		AssignACMG constructor

		Parameters
		----------
		variant_info: VariantInfo
			object containing variant basic info
		pathogenicity_odds_very_strong : float
			odds of pathogenicity for very strong evidence (OPVSt in Tavtigian et al. 2018)
		scaling_factor : float
			scaling factor to connect pathogenicity odds of different strength levels
		pathogenicity_prior : float
			pathogenicity prior
		"""
		self.logger = logging.getLogger("GenOtoScope_Classify.AssignACMG")
		self.logger.debug("Initialize class to assign ACMG class based on met criteria")
		self.variant_info = variant_info
		self.pathogenicity_odds_very_strong = pathogenicity_odds_very_strong
		self.scaling_factor = scaling_factor
		self.pathogenicity_prior = pathogenicity_prior

		assignments_preferred_inher = {"pathogenic_very_strong": {"PVS1": 0},
		                               "pathogenic_strong": {"PVS1_Strong": 0, "PS1": 0, "PM5_Strong": 0},
		                               "pathogenic_moderate": {"PVS1_Moderate": 0, "PM1": 0, "PM2": 0, "PM4": 0,
		                                                       "PM5_Moderate": 0},
		                               "pathogenic_supporting": {"PVS1_Supporting": 0, "PM2_Supporting": 0, "PP3": 0},
		                               "benign_stand_alone": {"BA1": 0}, "benign_strong": {"BS1": 0},
		                               "benign_supporting": {"BS1_Supporting": 0, "BP3": 0, "BP4": 0, "BP7": 0}}
		assignments_alternative_inher = {"pathogenic_very_strong": {"PVS1": 0},
		                                 "pathogenic_strong": {"PVS1_Strong": 0, "PS1": 0, "PM5_Strong": 0},
		                                 "pathogenic_moderate": {"PVS1_Moderate": 0, "PM1": 0, "PM2": 0, "PM4": 0,
		                                                         "PM5_Moderate": 0},
		                                 "pathogenic_supporting": {"PVS1_Supporting": 0, "PM2_Supporting": 0, "PP3": 0},
		                                 "benign_stand_alone": {"BA1": 0}, "benign_strong": {"BS1": 0},
		                                 "benign_supporting": {"BS1_Supporting": 0, "BP3": 0, "BP4": 0, "BP7": 0}}
		self.assignments = {"preferred": assignments_preferred_inher, "alternative": assignments_alternative_inher}

	@staticmethod
	def convert_assignment2int(rule_assignment):
		"""
		Convert rule assignment to integer assignment

		Parameters
		----------
		rule_assignment : str

		Returns
		-------
		int
			integer assignment of rule
		"""
		if rule_assignment == "True":
			return 1
		elif "PVS1" in rule_assignment:
			return 1
		elif "PM5" in rule_assignment:
			return 1
		elif rule_assignment == "False":
			return 0
		elif rule_assignment == "NA":
			return 0
		else:
			return -1

	def parse_assigned_rules(self, assigned_rules):
		"""
		Parse assigned rules to save met rules

		Parameters
		----------
		assigned_rules : str
			assignment of all rules for current variant

		Returns
		-------
		dict of str to int
			dictionary of parsed rules for preferred inheritance, key: rule name, value: rule assignment as integer
		str
			preferred inheritance
		dict of str to int
			dictionary of parsed rules for alternative inheritance, key: rule name, value: rule assignment as integer
		str
			alternative inheritance
		"""
		self.logger.debug("Parse assigned rules")
		preferred_inheritance, alternative_inheritance = None, None
		rules_preferred, rules_alternative = {}, {}
		is_preferred_seen = False
		rules = assigned_rules.split("||")
		for rule in rules:
			self.logger.debug("Split examined rule: {}".format(rule))
			if "NGSD" in rule or "AF" in rule or "QUAL" in rule:
				# exclude not ACMG rules from parsed output
				continue
			elif "=" in rule:
				### ### ###
				# extract inheritance related rules
				### ### ###
				if not is_preferred_seen:
					preferred_inheritance, rules_inheritance = rule.split("=")
					for rule_inheritance in rules_inheritance.split("|"):
						rule_name, rule_assignment = split_name_assignment(rule_inheritance)
						rules_preferred[rule_name] = AssignACMG.convert_assignment2int(rule_assignment)
						self.logger.debug(
							"parsed rule: name: {}, assignment: {}".format(rule_name, rules_preferred[rule_name]))
					is_preferred_seen = True
				else:
					alternative_inheritance, rules_inheritance = rule.split("=")
					for rule_inheritance in rules_inheritance.split("|"):
						rule_name, rule_assignment = split_name_assignment(rule_inheritance)
						rules_alternative[rule_name] = AssignACMG.convert_assignment2int(rule_assignment)
						self.logger.debug(
							"parsed rule: name:{}, assignment: {}".format(rule_name, rules_alternative[rule_name]))
			else:
				rule_name, rule_assignment = split_name_assignment(rule)
				rules_preferred[rule_name] = AssignACMG.convert_assignment2int(rule_assignment)
				self.logger.debug(
					"parsed rule: name: {}, assignment: {}".format(rule_name, rules_preferred[rule_name]))
				if alternative_inheritance:
					rules_alternative[rule_name] = AssignACMG.convert_assignment2int(rule_assignment)
		# assert that preferred inheritance is always parsed
		try:
			assert preferred_inheritance is not None
		except AssertionError:
			self.logger.error(
				"Parsing ACMG rules resulted to no preferred inheritance\n=> reference pos: {}".format(
					self.variant_info.to_string()), exc_info=True)

		self.logger.debug("Parsed rules preferred inheritance: {}".format(rules_preferred))
		self.logger.debug("Parsed rules alternative inheritance: {}".format(rules_alternative))

		return rules_preferred, preferred_inheritance, rules_alternative, alternative_inheritance

	def print_current_rules_assignments(self, selected_inheritance):
		"""
		Print current assignments of rules

		Parameters
		----------
		selected_inheritance : str
			print rule assignments for selected inheritance mode (preferred or alternative)

		Returns
		-------
		None
		"""
		self.logger.debug("### Current rules assignments ###")
		self.logger.debug("### Inheritance= {}           ###".format(selected_inheritance))
		self.logger.debug("###        Pathogenic         ###")
		self.logger.debug("Very strong: {}".format(self.assignments[selected_inheritance]["pathogenic_very_strong"]))
		self.logger.debug("Strong: {}".format(self.assignments[selected_inheritance]["pathogenic_strong"]))
		self.logger.debug("Moderate: {}".format(self.assignments[selected_inheritance]["pathogenic_moderate"]))
		self.logger.debug("Supporting: {}".format(self.assignments[selected_inheritance]["pathogenic_supporting"]))
		self.logger.debug("###          Benign          ###")
		self.logger.debug("Stand-alone: {}".format(self.assignments[selected_inheritance]["benign_stand_alone"]))
		self.logger.debug("Strong: {}".format(self.assignments[selected_inheritance]["benign_strong"]))
		self.logger.debug("Supporting: {}".format(self.assignments[selected_inheritance]["benign_supporting"]))
		self.logger.debug("###          #####           ###")

	def update_triggered_rules_counts(self, rules, selected_inheritance):
		"""
		Update counts of triggered ACMG rules for selected inheritance mode

		Parameters
		----------
		rules : dict of str : str
			dictionary of parsed ACMG rules, key: rule name, value: assignment result
		selected_inheritance : str
			selected inheritance mode to update rules count for (preferred, alternative)

		Returns
		-------
		None
		"""

		self.logger.debug("Update triggered rules for inheritance= {}".format(selected_inheritance))
		for rule_name, rule_result in rules.items():
			if rule_name in self.assignments[selected_inheritance]["pathogenic_very_strong"].keys():
				self.assignments[selected_inheritance]["pathogenic_very_strong"][rule_name] = rule_result
			elif rule_name in self.assignments[selected_inheritance]["pathogenic_strong"].keys():
				self.assignments[selected_inheritance]["pathogenic_strong"][rule_name] = rule_result
			elif rule_name in self.assignments[selected_inheritance]["pathogenic_moderate"].keys():
				self.assignments[selected_inheritance]["pathogenic_moderate"][rule_name] = rule_result
			elif rule_name in self.assignments[selected_inheritance]["pathogenic_supporting"].keys():
				self.assignments[selected_inheritance]["pathogenic_supporting"][rule_name] = rule_result
			elif rule_name in self.assignments[selected_inheritance]["benign_stand_alone"].keys():
				self.assignments[selected_inheritance]["benign_stand_alone"][rule_name] = rule_result
			elif rule_name in self.assignments[selected_inheritance]["benign_strong"].keys():
				self.assignments[selected_inheritance]["benign_strong"][rule_name] = rule_result
			elif rule_name in self.assignments[selected_inheritance]["benign_supporting"].keys():
				self.assignments[selected_inheritance]["benign_supporting"][rule_name] = rule_result
		self.print_current_rules_assignments(selected_inheritance)

	def deactivate_insignificant_pathogenic_rules(self, rules_comments, selected_inheritance):
		"""
		Deactivate automatically triggered pathogenic rules, that may be insignificant,
		aiming to enable benign classification with higher frequency

		Parameters
		----------
		rules_comments : str
			genotoscope comments for triggering rules
		selected_inheritance : str
			inheritance to deactivate pathogenic rules for

		Returns
		-------
		None
		"""
		self.logger.debug("Deactivate insignificant pathogenic rules, to enable benign classification")
		### ### ###
		# Deactivate pathogenic rules in favour of (potential) benign classification
		# 1) deactivate PM2 or PM2_Supporting triggered rule if variant not in gnomAD and not other pathogenic rules were triggered
		# 2) deactivate PP3 if it is the only pathogenic rule that was activated
		### ### ###

		# 1) deactivate PM2 or PM2_Supporting
		var_in_gnomad = True
		if "Inheritance_specific_rules: Not present in BOTH gnomAD subpopulations and general gnomAD, assume MAF=0" in str(
				rules_comments):
			var_in_gnomad = False

		if not var_in_gnomad:
			if sum(self.assignments[selected_inheritance]["pathogenic_very_strong"].values()) == 0 and sum(
					self.assignments[selected_inheritance]["pathogenic_strong"].values()) == 0 and (
					sum(self.assignments[selected_inheritance]["pathogenic_moderate"].values()) + sum(
					self.assignments[selected_inheritance]["pathogenic_supporting"].values()) == 1):
				if self.assignments[selected_inheritance]["pathogenic_moderate"]["PM2"] == 1:
					self.logger.debug("Variant not in gnomAD, deactivate PM2 rule")
					self.assignments[selected_inheritance]["pathogenic_moderate"]["PM2"] = 0
				if self.assignments[selected_inheritance]["pathogenic_supporting"]["PM2_Supporting"] == 1:
					self.logger.debug("Variant not in gnomAD, deactivate PM2_Supporting rule")
					self.assignments[selected_inheritance]["pathogenic_supporting"]["PM2_Supporting"] = 0

	def compute_probability(self, selected_inheritance, inheritance_value):
		"""
		Compute posterior pathogenicity probability for current variant
		based on equations 4 and 5 of Tavtigian et al. 2018

		Parameters
		----------
		gene_preferred_inheritance_mode: str
			preferred inheritance mode of variant-affected gene
		gene_alternative_inheritance_mode: str
			alternative inheritance mode of variant-affected gene

		Returns
		-------
		float
			Pathogenicity posterior probability
		"""
		self.logger.debug("Compute posterior probability of pathogenicity")
		if sum(self.assignments[selected_inheritance]["benign_stand_alone"].values()) == 1:
			return 0
		elif inheritance_value == "AR" and (
				sum(self.assignments[selected_inheritance]["pathogenic_very_strong"].values()) == 1 and
				self.assignments[selected_inheritance]["pathogenic_supporting"]["PM2_Supporting"] == 1):
			# first new combination of criteria for hearing loss
			# Please see "Results - Summary of specification" of Oza et al. (2018)
			evidence_exponent = sum(self.assignments[selected_inheritance]["pathogenic_very_strong"].values()) + (
					self.assignments[selected_inheritance]["pathogenic_supporting"][
						"PM2_Supporting"] / self.scaling_factor ** 3)
		elif sum(self.assignments[selected_inheritance]["benign_strong"].values()) == 1:
			# second new combination of criteria for hearing loss
			# Please see "Results - Summary of specification" of Oza et al. (2018)
			evidence_exponent = - (
					sum(self.assignments[selected_inheritance]["benign_strong"].values()) / self.scaling_factor)
		else:
			evidence_exponent = (sum(self.assignments[selected_inheritance][
				                         "pathogenic_supporting"].values()) / self.scaling_factor ** 3) + (
					                    sum(self.assignments[selected_inheritance][
						                        "pathogenic_moderate"].values()) / self.scaling_factor ** 2) + (
					                    sum(self.assignments[selected_inheritance][
						                        "pathogenic_strong"].values()) / self.scaling_factor) + sum(
				self.assignments[selected_inheritance]["pathogenic_very_strong"].values()) - (
					                    sum(self.assignments[selected_inheritance][
						                        "benign_supporting"].values()) / self.scaling_factor ** 3) - (
					                    sum(self.assignments[selected_inheritance][
						                        "benign_strong"].values()) / self.scaling_factor)
		pathogenicity_odds = self.pathogenicity_odds_very_strong ** evidence_exponent
		posterior_pathogenicity_prob = (pathogenicity_odds * self.pathogenicity_prior) / (
				(pathogenicity_odds - 1) * self.pathogenicity_prior + 1)
		try:
			assert 0 <= posterior_pathogenicity_prob <= 1
		except AssertionError:
			self.logger.error(
				"Posterior pathogenicity probability should be in the [0,1] range\n=> reference pos: {}".format(
					self.variant_info.to_string()), exc_info=True)
		return posterior_pathogenicity_prob

	def prepare4classify(self, examined_rules, rules_comments):
		"""
		Preprocess rules assignments to prepare for classification

		Parameters
		----------
		examined_rules : str
			examined ACMG rules for variant row
		rules_comments: str
			comments for triggering ACMG rules for variant row

		Returns
		-------
		str
			preferred inheritance mode
		str
			alternative inheritance mode
		"""
		self.logger.debug("Prepare for classification by parsing rules assignments")
		### ### ### ###
		# 1. parse rules for available inheritance modes
		# 2. update AssignACMG counts for all parsed rules of each inheritance
		# 3. deactivate not significant pathogenic rules
		### ### ### ###

		# 1) parse column with triggered rules
		parsed_rules_preferred, preferred_inheritance, parsed_rules_alternative, alternative_inheritance = self.parse_assigned_rules(
			examined_rules)

		# 2) update class counts for all parsed rules (per inheritance)
		self.update_triggered_rules_counts(parsed_rules_preferred, "preferred")
		if alternative_inheritance:
			self.update_triggered_rules_counts(parsed_rules_alternative, "alternative")

		# 3) deactivate not significant pathogenic rules
		self.deactivate_insignificant_pathogenic_rules(rules_comments, "preferred")
		if alternative_inheritance:
			self.deactivate_insignificant_pathogenic_rules(rules_comments, "alternative")

		return preferred_inheritance, alternative_inheritance

	def compute_categorical(self, selected_inheritance, inheritance_value):
		"""
		Compute categorical ACMG class (pathogenic, likely pathogenic, benign, likely benign and uncertain significance)
		based on the documented recommendation on how to combine triggered rules by Oza et al. 2018

		Parameters
		----------
		selected_inheritance : str
			selected inheritance (preferred or alternative)
		inheritance_value : str
			value of selected inheritance mode

		Returns
		-------
		str
			ACMG pathogenic class
		int
			ACMG pathogenic class (numeric)
		"""
		self.logger.debug("Compute ACMG pathogenicity class")
		self.logger.debug("Gene's {} inheritance: {}".format(selected_inheritance, inheritance_value))
		if sum(self.assignments[selected_inheritance]["benign_stand_alone"].values()) == 1 or sum(
				self.assignments[selected_inheritance]["benign_strong"].values()) >= 2:
			return "benign", 1
		elif (sum(self.assignments[selected_inheritance]["pathogenic_very_strong"].values()) >= 1 or sum(
				self.assignments[selected_inheritance]["pathogenic_strong"].values()) >= 1 or sum(
			self.assignments[selected_inheritance]["pathogenic_moderate"].values()) >= 1 or sum(
			self.assignments[selected_inheritance]["pathogenic_supporting"].values()) >= 1) and (
				sum(self.assignments[selected_inheritance]["benign_strong"].values()) >= 1 or sum(
			self.assignments[selected_inheritance]["benign_supporting"].values()) >= 1):
			# both pathogenic and benign criteria rules are triggered
			self.logger.debug("both pathogenic and benign criteria are triggered")
			return "uncertain significance", 3
		elif sum(self.assignments[selected_inheritance]["pathogenic_very_strong"].values()) == 1 and (
				sum(self.assignments[selected_inheritance]["pathogenic_strong"].values()) >= 1 or sum(
			self.assignments[selected_inheritance]["pathogenic_moderate"].values()) >= 2 or (
						sum(self.assignments[selected_inheritance]["pathogenic_moderate"].values()) == 1 and sum(
					self.assignments[selected_inheritance]["pathogenic_supporting"].values()) == 1) or sum(
			self.assignments[selected_inheritance]["pathogenic_supporting"].values()) >= 2):
			return "pathogenic", 5
		elif sum(self.assignments[selected_inheritance]["pathogenic_strong"].values()) >= 2:
			return "pathogenic", 5
		elif sum(self.assignments[selected_inheritance]["pathogenic_strong"].values()) == 1 and (
				sum(self.assignments[selected_inheritance]["pathogenic_moderate"].values()) >= 3 or (
				sum(self.assignments[selected_inheritance]["pathogenic_moderate"].values()) == 2 and sum(
			self.assignments[selected_inheritance]["pathogenic_supporting"].values()) >= 2) or (
						sum(self.assignments[selected_inheritance]["pathogenic_moderate"].values()) == 1 and sum(
					self.assignments[selected_inheritance]["pathogenic_supporting"].values()) >= 4)):
			return "pathogenic", 5
		elif inheritance_value == "AR" and (
				sum(self.assignments[selected_inheritance]["pathogenic_very_strong"].values()) == 1 and
				self.assignments[selected_inheritance]["pathogenic_supporting"]["PM2_Supporting"] == 1):
			# first new combination of criteria for hearing loss
			# Please see "Results - Summary of specification" of Oza et al. (2018)
			return "likely pathogenic", 4
		elif sum(self.assignments[selected_inheritance]["pathogenic_very_strong"].values()) == 1 and sum(
				self.assignments[selected_inheritance]["pathogenic_moderate"].values()) == 1:
			return "likely pathogenic", 4
		elif sum(self.assignments[selected_inheritance]["pathogenic_strong"].values()) == 1 and 1 <= sum(
				self.assignments[selected_inheritance]["pathogenic_moderate"].values()) <= 2:
			return "likely pathogenic", 4
		elif sum(self.assignments[selected_inheritance]["pathogenic_strong"].values()) == 1 and sum(
				self.assignments[selected_inheritance]["pathogenic_supporting"].values()) >= 2:
			return "likely pathogenic", 4
		elif sum(self.assignments[selected_inheritance]["pathogenic_moderate"].values()) >= 3:
			return "likely pathogenic", 4
		elif sum(self.assignments[selected_inheritance]["pathogenic_moderate"].values()) == 2 and sum(
				self.assignments[selected_inheritance]["pathogenic_supporting"].values()) >= 2:
			return "likely pathogenic", 4
		elif sum(self.assignments[selected_inheritance]["pathogenic_moderate"].values()) == 1 and sum(
				self.assignments[selected_inheritance]["pathogenic_supporting"].values()) >= 4:
			return "likely pathogenic", 4
		elif (sum(self.assignments[selected_inheritance]["benign_strong"].values()) == 1 and sum(
				self.assignments[selected_inheritance]["benign_supporting"].values()) == 1) or sum(
			self.assignments[selected_inheritance]["benign_supporting"].values()) >= 2:
			return "likely benign", 2
		elif sum(self.assignments[selected_inheritance]["benign_strong"].values()) == 1:
			# second new combination of criteria for hearing loss
			# Please see "Results - Summary of specification" of Oza et al. (2018)
			return "likely benign", 2
		else:
			# if none of the combinations is met then set "uncertain significance" class
			return "uncertain significance", 3
