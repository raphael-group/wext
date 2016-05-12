subroutine seed_prng(n)

    implicit none

    integer, intent(in) :: n
    integer :: m
    integer, allocatable :: x(:)

    call random_seed(size = m)
    allocate(x(m))
    x = n
    call random_seed(put = x)
    deallocate(x)

end subroutine seed_prng

subroutine bipartite_edge_swap(permuted_edge_list, edge_list, max_swaps, max_tries, seed, verbose, m, n, num_edges)

    implicit none

    integer, intent(in) :: max_swaps, max_tries, seed, verbose, m, n, num_edges
    integer, intent(in) :: edge_list(num_edges, 2)
    integer, intent(out) :: permuted_edge_list(num_edges, 2)

    logical :: permuted_adjacency_matrix(m, n)
    integer :: i, j, k, l, q, r, s, t, swaps, tries, modulo_tries, max_modulo_tries
    double precision :: x(max(min(max_swaps/100, 2**16), 1), 2)

    ! Initialize variables.
    permuted_edge_list = edge_list
    permuted_adjacency_matrix = .false.

    do k=1,num_edges
        i = edge_list(k, 1)
        j = edge_list(k, 2)
        permuted_adjacency_matrix(i, j) = .true.
    end do

    swaps = 0
    tries = 0

    ! Generate more random numbers at a time for speed.
    modulo_tries = max(min(max_swaps/100, 2**14), 1)
    max_modulo_tries = modulo_tries+1

    ! Seed random number generator.
    call seed_prng(seed)

    ! Repeat until reaching number of swaps or tries.
    do while ((swaps<max_swaps) .and. (tries<max_tries))

        ! Increment tries.
        tries = tries + 1
        modulo_tries = modulo_tries + 1

        ! Find two edges uniformly at random.
        if (modulo_tries==max_modulo_tries) then
            call random_number(x)
            modulo_tries = 1
        end if

        i = int(num_edges*x(modulo_tries, 1))+1
        j = int(num_edges*x(modulo_tries, 2))+1

        ! Try again if edges are same edge.
        if (i==j) then
            cycle
        end if

        ! Extract vertices from edges.
        q = permuted_edge_list(i, 1)
        r = permuted_edge_list(i, 2)
        s = permuted_edge_list(j, 1)
        t = permuted_edge_list(j, 2)

        ! Try again if swapped vertices already connected.
        if (permuted_adjacency_matrix(q, t) .or. permuted_adjacency_matrix(s, r)) then
            cycle
        end if

        ! Perform double edge swap: (q, r) and (s, t) to (q, t) and (s, r).
        permuted_edge_list(i, 2) = t
        permuted_edge_list(j, 2) = r

        permuted_adjacency_matrix(q, r) = .false.
        permuted_adjacency_matrix(s, t) = .false.
        permuted_adjacency_matrix(q, t) = .true.
        permuted_adjacency_matrix(s, r) = .true.

        ! Increment swaps.
        swaps = swaps + 1

    end do

    ! If verbose, then show information about swaps.
    if (verbose==1) then
        l = 0
        do k=1,num_edges
            i = edge_list(k, 1)
            j = edge_list(k, 2)
            if (permuted_adjacency_matrix(i, j)) then
                l = l+1
            end if
        end do

        print *, "Number of edges:          ", num_edges
        print *, "Number of preserved edges:", l
        print *, "Number of swaps:          ", swaps
        print *, "Maximum number of swaps:  ", max_swaps
        print *, "Number of tries:          ", tries
        print *, "Maximum number of tries:  ", max_tries

    end if

end subroutine bipartite_edge_swap
