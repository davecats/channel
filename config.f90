module config
  use, intrinsic :: iso_c_binding, only: C_DOUBLE, C_INT
  implicit none

  type :: config_entry
    character(len=64) :: section = ""
    character(len=64) :: key = ""
    character(len=256) :: value = ""
  end type config_entry

  type :: ini_config
    type(config_entry), allocatable :: entries(:)
  end type ini_config

contains

  subroutine read_ini_file(filename, cfg)
    character(len=*), intent(in) :: filename
    type(ini_config), intent(out) :: cfg

    integer :: unit, stat
    logical :: exists
    character(len=512) :: line
    character(len=64) :: section, key
    character(len=256) :: value

    if (allocated(cfg%entries)) deallocate (cfg%entries)

    inquire (file=trim(filename), exist=exists)
    if (.not. exists) then
      print *, "error: input file not found:", trim(filename)
      error stop "missing input file"
    end if

    open (newunit=unit, file=trim(filename), status="old", action="read", iostat=stat)
    if (stat /= 0) then
      print *, "error: could not open input file:", trim(filename)
      error stop "could not open input file"
    end if

    section = ""
    do
      read (unit, '(A)', iostat=stat) line
      if (stat /= 0) exit

      call strip_comment(line)
      line = adjustl(line)
      if (len_trim(line) == 0) cycle

      if (line(1:1) == "[") then
        call parse_section(line, section)
        cycle
      end if

      call split_key_value(line, key, value)
      if (len_trim(key) == 0) cycle
      call append_entry(cfg, lower(trim(section)), lower(trim(key)), trim(value))
    end do

    close (unit)
  end subroutine read_ini_file

  subroutine require_integer(cfg, section, key, target)
    type(ini_config), intent(in) :: cfg
    character(len=*), intent(in) :: section, key
    integer(C_INT), intent(out) :: target

    logical :: found

    call get_integer(cfg, section, key, target, found)
    if (.not. found) call missing_key(section, key)
  end subroutine require_integer

  subroutine require_real(cfg, section, key, target)
    type(ini_config), intent(in) :: cfg
    character(len=*), intent(in) :: section, key
    real(C_DOUBLE), intent(out) :: target

    logical :: found

    call get_real(cfg, section, key, target, found)
    if (.not. found) call missing_key(section, key)
  end subroutine require_real

  subroutine require_logical(cfg, section, key, target)
    type(ini_config), intent(in) :: cfg
    character(len=*), intent(in) :: section, key
    logical, intent(out) :: target

    logical :: found

    call get_logical(cfg, section, key, target, found)
    if (.not. found) call missing_key(section, key)
  end subroutine require_logical

  subroutine require_real_vector(cfg, section, key, target)
    type(ini_config), intent(in) :: cfg
    character(len=*), intent(in) :: section, key
    real(C_DOUBLE), intent(out) :: target(:)

    logical :: found

    call get_real_vector(cfg, section, key, target, found)
    if (.not. found) call missing_key(section, key)
  end subroutine require_real_vector

  subroutine get_integer(cfg, section, key, target, found)
    type(ini_config), intent(in) :: cfg
    character(len=*), intent(in) :: section, key
    integer(C_INT), intent(inout) :: target
    logical, intent(out) :: found

    integer :: idx, stat, parsed
    character(len=256) :: value

    call lookup_value(cfg, section, key, value, found, idx)
    if (.not. found) return

    read (value, *, iostat=stat) parsed
    if (stat /= 0) then
      print *, "error: could not parse integer value for", trim(section)//"."//trim(key)
      error stop "invalid integer value"
    end if
    target = int(parsed, C_INT)
  end subroutine get_integer

  subroutine get_real(cfg, section, key, target, found)
    type(ini_config), intent(in) :: cfg
    character(len=*), intent(in) :: section, key
    real(C_DOUBLE), intent(inout) :: target
    logical, intent(out) :: found

    integer :: idx, stat
    character(len=256) :: value

    call lookup_value(cfg, section, key, value, found, idx)
    if (.not. found) return

    read (value, *, iostat=stat) target
    if (stat /= 0) then
      print *, "error: could not parse real value for", trim(section)//"."//trim(key)
      error stop "invalid real value"
    end if
  end subroutine get_real

  subroutine get_logical(cfg, section, key, target, found)
    type(ini_config), intent(in) :: cfg
    character(len=*), intent(in) :: section, key
    logical, intent(inout) :: target
    logical, intent(out) :: found

    integer :: idx
    character(len=256) :: value

    call lookup_value(cfg, section, key, value, found, idx)
    if (.not. found) return

    call parse_logical_value(value, target)
  end subroutine get_logical

  subroutine get_string(cfg, section, key, target, found)
    type(ini_config), intent(in) :: cfg
    character(len=*), intent(in) :: section, key
    character(len=*), intent(inout) :: target
    logical, intent(out) :: found

    integer :: idx
    character(len=256) :: value

    call lookup_value(cfg, section, key, value, found, idx)
    if (.not. found) return

    call clean_string(value, target)
  end subroutine get_string

  subroutine get_real_vector(cfg, section, key, target, found)
    type(ini_config), intent(in) :: cfg
    character(len=*), intent(in) :: section, key
    real(C_DOUBLE), intent(inout) :: target(:)
    logical, intent(out) :: found

    integer :: idx, stat
    character(len=256) :: value

    call lookup_value(cfg, section, key, value, found, idx)
    if (.not. found) return

    if (size(target) == 0) return

    read (value, *, iostat=stat) target
    if (stat /= 0) then
      print *, "error: could not parse real vector value for", trim(section)//"."//trim(key)
      error stop "invalid real vector value"
    end if
  end subroutine get_real_vector

  logical function has_section(cfg, section)
    type(ini_config), intent(in) :: cfg
    character(len=*), intent(in) :: section

    integer :: i
    character(len=64) :: section_l

    has_section = .false.
    if (.not. allocated(cfg%entries)) return

    section_l = lower(trim(section))
    do i = 1, size(cfg%entries)
      if (trim(cfg%entries(i)%section) == trim(section_l)) then
        has_section = .true.
        return
      end if
    end do
  end function has_section

  subroutine append_entry(cfg, section, key, value)
    type(ini_config), intent(inout) :: cfg
    character(len=*), intent(in) :: section, key, value

    type(config_entry), allocatable :: tmp(:)
    integer :: n

    if (.not. allocated(cfg%entries)) then
      allocate (cfg%entries(1))
      n = 0
    else
      n = size(cfg%entries)
      allocate (tmp(n))
      tmp = cfg%entries
      deallocate (cfg%entries)
      allocate (cfg%entries(n + 1))
      cfg%entries(1:n) = tmp
      deallocate (tmp)
    end if

    cfg%entries(n + 1)%section = ""
    cfg%entries(n + 1)%key = ""
    cfg%entries(n + 1)%value = ""
    cfg%entries(n + 1)%section = trim(section)
    cfg%entries(n + 1)%key = trim(key)
    cfg%entries(n + 1)%value = trim(value)
  end subroutine append_entry

  subroutine lookup_value(cfg, section, key, value, found, idx)
    type(ini_config), intent(in) :: cfg
    character(len=*), intent(in) :: section, key
    character(len=*), intent(out) :: value
    logical, intent(out) :: found
    integer, intent(out) :: idx

    character(len=64) :: section_l, key_l

    value = ""
    found = .false.
    idx = 0
    if (.not. allocated(cfg%entries)) return

    section_l = lower(trim(section))
    key_l = lower(trim(key))

    do idx = size(cfg%entries), 1, -1
      if (trim(cfg%entries(idx)%section) == trim(section_l) .and. &
          trim(cfg%entries(idx)%key) == trim(key_l)) then
        value = trim(cfg%entries(idx)%value)
        found = .true.
        return
      end if
    end do
  end subroutine lookup_value

  subroutine parse_section(line, section)
    character(len=*), intent(in) :: line
    character(len=*), intent(out) :: section
    integer :: last

    section = ""
    last = index(line, "]")
    if (last > 2) section = lower(trim(line(2:last - 1)))
  end subroutine parse_section

  subroutine split_key_value(line, key, value)
    character(len=*), intent(in) :: line
    character(len=*), intent(out) :: key, value
    integer :: eq

    key = ""
    value = ""
    eq = index(line, "=")
    if (eq <= 0) return

    key = trim(adjustl(line(:eq - 1)))
    value = trim(adjustl(line(eq + 1:)))
  end subroutine split_key_value

  subroutine strip_comment(line)
    character(len=*), intent(inout) :: line
    integer :: bang, hash, semicolon, cut

    bang = index(line, "!")
    hash = index(line, "#")
    semicolon = index(line, ";")
    cut = 0

    if (bang > 0) cut = bang
    if (hash > 0 .and. (cut == 0 .or. hash < cut)) cut = hash
    if (semicolon > 0 .and. (cut == 0 .or. semicolon < cut)) cut = semicolon
    if (cut > 0) line(cut:) = ""
  end subroutine strip_comment

  subroutine parse_logical_value(text, value)
    character(len=*), intent(in) :: text
    logical, intent(out) :: value

    character(len=256) :: cleaned

    call clean_string(text, cleaned)
    cleaned = lower(trim(cleaned))

    select case (trim(cleaned))
    case ("true", ".true.", "1", "yes", "on")
      value = .true.
    case ("false", ".false.", "0", "no", "off")
      value = .false.
    case default
      print *, "error: could not parse logical value:", trim(text)
      error stop "invalid logical value"
    end select
  end subroutine parse_logical_value

  subroutine clean_string(text, out)
    character(len=*), intent(in) :: text
    character(len=*), intent(out) :: out

    character(len=256) :: tmp
    integer :: n

    out = ""
    tmp = trim(adjustl(text))
    n = len_trim(tmp)
    if (n == 0) return

    if (n >= 2) then
      if ((tmp(1:1) == '"' .and. tmp(n:n) == '"') .or. &
          (tmp(1:1) == "'" .and. tmp(n:n) == "'")) then
        if (n > 2) out = tmp(2:n - 1)
        return
      end if
    end if

    out = tmp(:min(len(out), n))
  end subroutine clean_string

  function lower(text) result(out)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: out

    integer :: i, c

    do i = 1, len(text)
      c = iachar(text(i:i))
      if (c >= iachar("A") .and. c <= iachar("Z")) then
        out(i:i) = achar(c + iachar("a") - iachar("A"))
      else
        out(i:i) = text(i:i)
      end if
    end do
  end function lower

  subroutine missing_key(section, key)
    character(len=*), intent(in) :: section, key

    print *, "error: missing required config key", trim(section)//"."//trim(key)
    error stop "missing config key"
  end subroutine missing_key

end module config
