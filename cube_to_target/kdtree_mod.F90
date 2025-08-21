module kdtree_mod
    implicit none
    
    type :: kdtree_node
        integer :: point_idx
        real(8) :: coords(2)  ! lon, lat
        type(kdtree_node), pointer :: left => null()
        type(kdtree_node), pointer :: right => null()
        integer :: split_dim  ! 1 for lon, 2 for lat
    end type kdtree_node
    
    type :: kdtree
        type(kdtree_node), pointer :: root => null()
        integer :: num_points
        real(8), allocatable :: points(:,:)  ! (2, num_points) - lon, lat
        integer, allocatable :: indices(:)   ! original indices
        logical, allocatable :: valid(:)     ! valid flags
    end type kdtree
    
contains

    ! Build k-d tree from valid points
    subroutine build_kdtree(tree, target_center_lon, target_center_lat, valid_cells)
        type(kdtree), intent(out) :: tree
        real(8), intent(in) :: target_center_lon(:), target_center_lat(:)
        logical, intent(in) :: valid_cells(:)
        
        integer :: i, valid_count, idx
        integer, allocatable :: temp_indices(:)
        
        ! Count valid points
        valid_count = count(valid_cells)
        tree%num_points = valid_count
        
        if (valid_count == 0) return
        
        ! Allocate arrays
        allocate(tree%points(2, valid_count))
        allocate(tree%indices(valid_count))
        allocate(tree%valid(valid_count))
        allocate(temp_indices(valid_count))
        
        ! Copy valid points
        idx = 1
        do i = 1, size(valid_cells)
            if (valid_cells(i)) then
                tree%points(1, idx) = target_center_lon(i)
                tree%points(2, idx) = target_center_lat(i)
                tree%indices(idx) = i
                tree%valid(idx) = .true.
                temp_indices(idx) = idx
                idx = idx + 1
            end if
        end do
        
        ! Build the tree
        tree%root => build_node(tree, temp_indices, 1, valid_count, 1)
        
        deallocate(temp_indices)
    end subroutine build_kdtree

    ! Recursive function to build k-d tree nodes
    recursive function build_node(tree, indices, start_idx, end_idx, depth) result(node)
        type(kdtree), intent(inout) :: tree
        integer, intent(inout) :: indices(:)
        integer, intent(in) :: start_idx, end_idx, depth
        type(kdtree_node), pointer :: node
        
        integer :: median_idx, split_dim
        
        if (start_idx > end_idx) then
            node => null()
            return
        end if
        
        allocate(node)
        
        ! Determine split dimension (alternate between lon and lat)
        split_dim = mod(depth - 1, 2) + 1
        node%split_dim = split_dim
        
        if (start_idx == end_idx) then
            ! Leaf node
            node%point_idx = indices(start_idx)
            node%coords(1) = tree%points(1, indices(start_idx))
            node%coords(2) = tree%points(2, indices(start_idx))
            node%left => null()
            node%right => null()
        else
            ! Find median and partition
            median_idx = (start_idx + end_idx) / 2
            call quickselect(tree, indices, start_idx, end_idx, median_idx, split_dim)
            
            node%point_idx = indices(median_idx)
            node%coords(1) = tree%points(1, indices(median_idx))
            node%coords(2) = tree%points(2, indices(median_idx))
            
            ! Recursively build left and right subtrees
            node%left => build_node(tree, indices, start_idx, median_idx - 1, depth + 1)
            node%right => build_node(tree, indices, median_idx + 1, end_idx, depth + 1)
        end if
    end function build_node

    ! Quickselect algorithm to find median
    recursive subroutine quickselect(tree, indices, left, right, k, dim)
        type(kdtree), intent(in) :: tree
        integer, intent(inout) :: indices(:)
        integer, intent(in) :: left, right, k, dim
        
        integer :: pivot_idx 
        
        if (left >= right) return
        
        pivot_idx = partition(tree, indices, left, right, dim)
        
        if (k == pivot_idx) then
            return
        else if (k < pivot_idx) then
            call quickselect(tree, indices, left, pivot_idx - 1, k, dim)
        else
            call quickselect(tree, indices, pivot_idx + 1, right, k, dim)
        end if
    end subroutine quickselect

    ! Partition function for quickselect
    function partition(tree, indices, left, right, dim) result(pivot_idx)
        type(kdtree), intent(in) :: tree
        integer, intent(inout) :: indices(:)
        integer, intent(in) :: left, right, dim
        integer :: pivot_idx
        
        real(8) :: pivot_value
        integer :: i, store_idx, temp
        
        ! Use rightmost element as pivot
        pivot_value = tree%points(dim, indices(right))
        store_idx = left
        
        do i = left, right - 1
            if (tree%points(dim, indices(i)) <= pivot_value) then
                ! Swap
                temp = indices(i)
                indices(i) = indices(store_idx)
                indices(store_idx) = temp
                store_idx = store_idx + 1
            end if
        end do
        
        ! Move pivot to final position
        temp = indices(right)
        indices(right) = indices(store_idx)
        indices(store_idx) = temp
        
        pivot_idx = store_idx
    end function partition

    ! Find nearest neighbor
    function find_nearest_neighbor_kdtree(tree, query_lon, query_lat, exclude_idx) result(nearest_idx)
        type(kdtree), intent(in) :: tree
        real(8), intent(in) :: query_lon, query_lat
        integer, intent(in), optional :: exclude_idx
        integer :: nearest_idx
        
        real(8) :: query_point(2)
        real(8) :: best_dist_sq
        integer :: best_idx, exclude_id
        
        if (.not. associated(tree%root)) then
            nearest_idx = -1
            return
        end if
        
        query_point(1) = query_lon
        query_point(2) = query_lat
        best_dist_sq = huge(1.0d0)
        best_idx = -1
        
        exclude_id = -1
        if (present(exclude_idx)) exclude_id = exclude_idx
        
        call search_nearest(tree, tree%root, query_point, best_dist_sq, best_idx, exclude_id)
        
        if (best_idx > 0) then
            nearest_idx = tree%indices(best_idx)
        else
            nearest_idx = -1
        end if
    end function find_nearest_neighbor_kdtree

    ! Recursive nearest neighbor search
    recursive subroutine search_nearest(tree, node, query_point, best_dist_sq, best_idx, exclude_idx)
        type(kdtree), intent(in) :: tree
        type(kdtree_node), pointer, intent(in) :: node
        real(8), intent(in) :: query_point(2)
        real(8), intent(inout) :: best_dist_sq
        integer, intent(inout) :: best_idx
        integer, intent(in) :: exclude_idx
        
        real(8) :: dist_sq, diff_sq
        type(kdtree_node), pointer :: first_child, second_child
        
        if (.not. associated(node)) return
        
        ! Calculate distance to current node
        dist_sq = (node%coords(1) - query_point(1))**2 + (node%coords(2) - query_point(2))**2
        
        ! Update best if this is better and not excluded
        if (dist_sq < best_dist_sq .and. tree%indices(node%point_idx) /= exclude_idx) then
            best_dist_sq = dist_sq
            best_idx = node%point_idx
        end if
        
        ! Determine which child to search first
        if (query_point(node%split_dim) <= node%coords(node%split_dim)) then
            first_child => node%left
            second_child => node%right
        else
            first_child => node%right
            second_child => node%left
        end if
        
        ! Search the closer child first
        call search_nearest(tree, first_child, query_point, best_dist_sq, best_idx, exclude_idx)
        
        ! Check if we need to search the other child
        diff_sq = (query_point(node%split_dim) - node%coords(node%split_dim))**2
        if (diff_sq < best_dist_sq) then
            call search_nearest(tree, second_child, query_point, best_dist_sq, best_idx, exclude_idx)
        end if
    end subroutine search_nearest

    ! Clean up k-d tree memory
    subroutine destroy_kdtree(tree)
        type(kdtree), intent(inout) :: tree
        
        if (associated(tree%root)) then
            call destroy_node(tree%root)
        end if
        
        if (allocated(tree%points)) deallocate(tree%points)
        if (allocated(tree%indices)) deallocate(tree%indices)
        if (allocated(tree%valid)) deallocate(tree%valid)
        
        tree%num_points = 0
    end subroutine destroy_kdtree

    ! Recursive node destruction
    recursive subroutine destroy_node(node)
        type(kdtree_node), pointer :: node
        
        if (.not. associated(node)) return
        
        call destroy_node(node%left)
        call destroy_node(node%right)
        deallocate(node)
        node => null()
    end subroutine destroy_node

end module kdtree_mod
