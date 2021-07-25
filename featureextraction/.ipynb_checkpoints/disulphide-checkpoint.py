def getPositionList(sequence):
    position_list = []

    for index, amino in enumerate(sequence):
        if amino == 'C':
            position_list.append(index+1)

    return position_list


def getDisulphideProperties(sequence):

    parallel_bond = 0
    alternate_bond = 0
    quad_core_bond = 0

    count_parallel_bond = 0
    count_alternate_bond = 0
    count_quad_core_bond = 0

    position_list = getPositionList(sequence)
    length = len(position_list)

    if length == 0 or length == 1:
        return (0, 0, 0)

    for index in range(1, length, 4):
        element1 = position_list[index]
        element2 = position_list[index-1]

        parallel_bond += (element1 - element2)
        count_parallel_bond += 1

        if index >= 2:
            element3 = position_list[index-2]

            if index != length-1:
                alternate_bond += (element1 - element3)
                count_alternate_bond += 1

            quad_core_bond += (element2 - element3)
            count_quad_core_bond += 1

        if index >= 3:
            element4 = position_list[index-3]

            if index != length-1:
                parallel_bond += (element3 - element4)
                quad_core_bond += (element1 - element4)

                count_quad_core_bond += 1
                count_parallel_bond += 1

            alternate_bond += (element2 - element4)
            count_alternate_bond += 1

    count_alternate_bond = 1 if count_alternate_bond == 0 else count_alternate_bond
    count_quad_core_bond = 1 if count_quad_core_bond == 0 else count_quad_core_bond

    return (
        parallel_bond / count_parallel_bond,
        alternate_bond / count_alternate_bond,
        quad_core_bond / count_quad_core_bond
    )
