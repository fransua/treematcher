def match(self_node, target_node):

	#status: self_node matches tearget_node
	status = seflf_node == target_node

	if self_node != target_node:
		if self_node.parent_node == '+' #means self_node can be skipped so treat it as a '+' too!
		status = True # because '+' alway matches


	# so 

	if status and self_node.has_children: 		 # the search can continue

		nodes_to_be_searched = [] 

		if target_node.children < self_node.children:

			if self_node == '+':  # self_node do not fit there, but it can matche later

				traverse_all_the_lower:

					if lower_node.has_the_right_number_of_children(): 		# right_number='>='
						
						lower_node.put_in_nodes_to_be_searched()

						if lower_node.has_sister():
							put_the_lower_node_SISTER_in_nodes_to_be_searched()		 # why? :/

							stop_traversing()

				if did_not_find_possible_matches():
					status = false

			if 	self_node != '+': 		# can not be skipped
				status = false


		if len(nodes_to_be_searched) == 0: 		# there is no starting point yet 
			self_node.put_in_nodes_to_be_searched()

		for  every nodes_to_be_searched:

			super_visor_count_vartiable = 0 		# Accumulator parameter

			if target_node.number_of_children >= self_node.number_of_children:

				super_visor_for_status = True

				for every_possible_compbination_of_children_sequence(): 		# guess is (A, B) or (B, A)					

					for every child:
						see_if_matches()

						# once leaves are met, lets "return" and see there is match. 


						if self_node_is_not_match() and self_node.can_be_skipped():
							ignore_self_node()

						else:
							super_visor_count_vartiable += 1 	# +1 for every checked node
							a_logical_mess_here()


					if super_visor_for_status = True and super_visor_count_vartiable > 0: 		#	has checked at least one node
						status = True 		# universal match found!!! 
						stop_searching()	#break
					else:
						continue_search() 	# status = False

				if lower_nodes_matched and super_visor_count_vartiable > 0: 		#	has checked at least one node
					stop_searching()	#break

	return status
	



def is_local_match(self_node, target_node):


	values_scope_is_"this"_plus_the_reference_to_target_node()

	if self_node.has_no_name() or self_node.can_be_skipped():
		no_constraint()

	if self_node.has_name():
		constraint = names_have_to_match()

	if self_node.name has('@'):
		constraint = replace(@, reference_to_node)


	try:
		if there_is_constraint():
			MATCH_FOUND = evaluate_constraint_with_names_taken_from_the_"this"_scope()
		else:
			MATCH_FOUND = True 

	except name or value error:
		can_not_be_fixed()

	except name_error:
		try:
			retry_using_roots_syntax()
			retrurn MATCH_FOUND
		except name_error:
			can_not_be_fixed()

	else:
		return MATCH_FOUND