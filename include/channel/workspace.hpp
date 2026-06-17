#pragma once

#include "channel/device_vector.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>

namespace channel {

enum class WorkspaceBackend {
  Native,
  UmpireSpace,
};

class WorkspaceArena;

class WorkspaceLease {
public:
  WorkspaceLease() = default;
  WorkspaceLease(WorkspaceArena* arena, std::string owner, std::size_t bytes, std::uint8_t* base);
  WorkspaceLease(const WorkspaceLease&) = delete;
  WorkspaceLease& operator=(const WorkspaceLease&) = delete;
  WorkspaceLease(WorkspaceLease&& other) noexcept;
  WorkspaceLease& operator=(WorkspaceLease&& other) noexcept;
  ~WorkspaceLease();

  [[nodiscard]] std::size_t bytes() const { return bytes_; }
  [[nodiscard]] std::size_t used() const { return cursor_; }
  [[nodiscard]] std::uint8_t* data() const { return base_; }

  template <class T>
  T* slice(std::size_t count) {
    const std::size_t aligned = align_offset(cursor_, alignof(T));
    const std::size_t nbytes = sizeof(T) * count;
    if (aligned + nbytes > bytes_) {
      throw std::runtime_error("WorkspaceLease slice exceeds lease capacity");
    }
    cursor_ = aligned + nbytes;
    return reinterpret_cast<T*>(base_ + aligned);
  }

private:
  friend class WorkspaceArena;
  static std::size_t align_offset(std::size_t offset, std::size_t alignment);
  void release();

  WorkspaceArena* arena_ = nullptr;
  std::string owner_;
  std::size_t bytes_ = 0;
  std::size_t cursor_ = 0;
  std::uint8_t* base_ = nullptr;
};

class WorkspaceArena {
public:
  explicit WorkspaceArena(WorkspaceBackend backend = WorkspaceBackend::UmpireSpace);
  WorkspaceArena(const WorkspaceArena&) = delete;
  WorkspaceArena& operator=(const WorkspaceArena&) = delete;
  WorkspaceArena(WorkspaceArena&&) = delete;
  WorkspaceArena& operator=(WorkspaceArena&&) = delete;
  ~WorkspaceArena();

  WorkspaceLease lease(std::size_t bytes, const std::string& owner);
  void reserve(std::size_t bytes);
  void release(const std::string& owner);

  [[nodiscard]] WorkspaceBackend backend() const { return backend_; }
  [[nodiscard]] std::size_t capacity_bytes() const { return capacity_; }
  [[nodiscard]] std::size_t high_water_bytes() const { return high_water_; }
  [[nodiscard]] bool owned() const { return owned_; }
  [[nodiscard]] const std::string& owner() const { return owner_; }

  static constexpr std::size_t alignment = 256;

private:
  [[nodiscard]] std::uint8_t* storage_data();

  WorkspaceBackend backend_;
  DeviceVector<std::uint8_t> storage_;
  std::size_t capacity_ = 0;
  std::size_t high_water_ = 0;
  bool owned_ = false;
  std::string owner_;
};

} // namespace channel
